# -*- coding: utf-8 -*-
# __author__ = 'xuting'
import web
import random
import json
import time
from ..core.basic import Basic
from mainapp.models.mongo.core.base import Base
from mainapp.libs.input_check import meta_check
from mainapp.models.mongo.meta import Meta
from mainapp.libs.signature import check_sig
from mainapp.models.mongo.distance_matrix import Distance
from mainapp.models.workflow import Workflow
from biocluster.core.function import filter_error_info
from biocluster.config import Config
# from mainapp.config.db import Config
import os


class MetaController(Base):

    def __init__(self, bind_object=None, instant=False):
        super(MetaController, self).__init__(bind_object)
        self._project_type = "meta"
        self._instant = instant
        self._post_data = None
        self._sheet_data = None
        self._return_msg = None
        # self.mongodb = Config().MONGODB
        self.meta = Meta()
        self.workflow_id = ""

    @property
    def data(self):
        """
        获取Post数据
        :return:
        """
        return self._post_data

    @property
    def return_msg(self):
        """
        获取Post数据

        :return:
        """
        return self._return_msg

    @property
    def instant(self):
        """
        任务是否是即时计算

        :return: bool
        """
        return self._instant

    @property
    def sheet_data(self):
        """
        获取运行流程所需的Json数据

        :return:
        """
        return self._sheet_data


    @check_sig
    @meta_check
    def POST(self):
        workflow_client = Basic(data=self.sheet_data, instant=self.instant)
        try:
            run_info = workflow_client.run()
            print("打印出run_info： {}".format(run_info))
            try:
                info_msg = json.loads(run_info['info'])  # info嵌套info的情况
                if info_msg['code'] == 'R001' and info_msg['info'] != "":
                    del info_msg['code']
                run_info.update(info_msg)
            except:
                run_info['info'] = run_info['info']  # info只包含字符的情况
            self._return_msg = workflow_client.return_msg
            print("修改后的run_info： {}".format(run_info))
            if not run_info['success']:
                self.roll_back()
                return {"success": False, "info": "Failed to submit your task because too many tasks are being "
                                                  "calculated on the cluster. Please try again!"}
            else:
                return run_info
        except Exception:
            self.roll_back()
            return {"success": False, "info": "Failed to submit your task because too many tasks are being "
                                              "calculated on the cluster. Please try again!"}

    def roll_back(self):
        """
        当任务投递失败时，如WPM服务出错时，主表写入start状态无法由API更新，此处进行更新

        :return:
        """
        print("INFO: 任务提交出错，尝试更新主表状态为failed。")
        try:
            update_info = json.loads(self.sheet_data['options']['update_info'])
            for i in update_info:
                if i == "batch_id":
                    continue
                self.meta.update_status_failed(update_info[i], i)
                print("INFO: 更新主表状态为failed成功: coll:{} _id:{}".format(update_info[i], i))
        except Exception as e:
            print('ERROR:尝试回滚主表状态为failed 失败:{}'.format(e))

    '''
    def get_main_info(self, main_id, collection_name):
        if isinstance(main_id, types.StringTypes):
            main_id = ObjectId(main_id)
        elif isinstance(main_id, ObjectId):
            main_id = main_id
        else:
            raise Exception("输入main_id参数必须为字符串或者ObjectId类型!")
        collection = self.db[collection_name]
        express_info = collection.find_one({'_id': main_id})
        return express_info
    '''

    def set_sheet_data(self, name, options, main_table_name, module_type="workflow", params=None, to_file=None,
                       main_id=None, collection_name=None):
        """
        设置运行所需的Json文档

        :param name: workflow/module/tool相对路径
        :param module_type: workflow/module/tool
        :param main_table_name: 交互分析项主表名称
        :param options: workflow/module/tool参数
        :param params: 交互分析主表params字段
        :param to_file: workflow/module/tool mongo数据转文件
        :param main_id:
        :param collection_name
        :return:
        """
        self._post_data = web.input()
        if not main_id:
            main_id = self.data.otu_id
            collection_name = 'sg_otu'
        table_info = Meta().get_main_info(main_id=main_id, collection_name=collection_name)
        # table_info = self.get_main_info(main_id=main_id, collection_name=collection_name)  # modified by guhaidong 20171026
        # print table_info
        project_sn = table_info["project_sn"]
        task_id = table_info["task_id"]
        new_task_id = self.get_new_id(task_id)
        self._sheet_data = {
            'id': new_task_id,
            "batch": False,
            'interaction': True,
            'name': name,  # 需要配置
            'type': module_type,  # 可以配置
            'client': self.data.client,
            'output': self._create_output_dir(task_id, main_table_name),
            'project_sn': project_sn,
            'IMPORT_REPORT_DATA': True,
            'UPDATE_STATUS_API': self._update_status_api(),
            'db_type': "meta",
            'options': options  # 需要配置
        }
        if "meta_base" in name:
            query_id = table_info["query_id"]
            self._sheet_data["query_id"] = query_id
        if self.instant:
            self._sheet_data["instant"] = True
        if params:
            self._sheet_data["params"] = params
        if to_file:
            self._sheet_data["to_file"] = to_file
        # print('Sheet_Data: {}'.format(self._sheet_data))
        # 此处不能print  sheet_data的量可能比较大，会造成uwsgi日志问题
        self.workflow_id = new_task_id
        self.meta_pipe()
        return self._sheet_data

    def meta_pipe(self):
        """
        一键化分析特殊处理
        """
        data = web.input()
        # print "data", data
        for i in ["batch_id"]:
            if not hasattr(data, i):
                return
        print "一键化投递任务{}: {}".format(i, getattr(data, i))
        if not hasattr(self.data, "batch_task_id"):
            print "NO BATCH_TASK_ID"
        else:
            print "BATCH_task_id " + self.data.batch_task_id
            self._sheet_data["batch_id"] = self.data.batch_task_id
        update_info = json.loads(self._sheet_data["options"]['update_info'])
        # update_info["meta_pipe_detail_id"] = data.meta_pipe_detail_id
        update_info["batch_id"] = data.batch_id
        self._sheet_data['options']["update_info"] = json.dumps(update_info)
        # if self._sheet_data['name'].strip().split(".")[-1] not in ["otu_subsample"]:
        #     self._instant = False
        #     self._sheet_data["instant"] = False
        self._instant = False
        self._sheet_data["instant"] = False
        self._sheet_data["batch"] = True

    def get_new_id(self, task_id, otu_id=None):
        """
        根据旧的ID生成新的workflowID，固定为旧的后面用“_”，添加两次随机数或者一次otu_id一次随机数
        """
        if otu_id:
            new_id = "{}_{}_{}".format(task_id, otu_id[-4:], random.randint(1, 10000))
        else:
            id_ = '%f' % time.time()
            ids = str(id_).strip().split(".")
            new_id = "{}_{}_{}".format(task_id, ids[0][5:], ids[1])  # 改成时间来命名workflow id
            # new_id = "{}_{}_{}".format(task_id, random.randint(1000, 10000), random.randint(1, 10000))
        workflow_module = Workflow()
        workflow_data = workflow_module.get_by_workflow_id(new_id)
        if len(workflow_data) > 0:
            return self.get_new_id(task_id, otu_id)
        return new_id

    def _create_output_dir(self, task_id, main_table_name):
        """
        根据主表名称，生成结果目录名称/上传路径

        modified by hongdongxuan 20170320
        """
        data = web.input()
        task_info = self.meta.get_task_info(task_id)
        if Config().RGW_ENABLE:
            part_dir = str(task_info['member_id']) + \
                          '/' + str(task_info['project_sn']) + '/' + \
                          task_id + '/interaction_results/' + main_table_name
            sanger_prefix = Config().get_project_region_bucket(project_type="meta")
            target_dir = os.path.join(sanger_prefix, 'files/' + part_dir + '/')  # 需要配置bucket, 改为对象存储路径 guhaidong 20180927
        else:
            client = data.client if hasattr(data, "client") else web.ctx.env.get('HTTP_CLIENT')
            if client == 'client01':
                target_dir = 'sanger'
            else:
                target_dir = 'tsanger'
            # target_dir += ':rerewrweset/files/' + str(task_info['member_id']) + \
            #               '/' + str(task_info['project_sn']) + '/' + \
            #               task_id + '/report_results/' + main_table_name
            target_dir += ':rerewrweset/files/' + str(task_info['member_id']) + \
                          '/' + str(task_info['project_sn']) + '/' + \
                          task_id + '/interaction_results/' + main_table_name  # zengjing 20170929 修改页面上meta多样性交互分析的结果文件夹名称为interaction_results
        return target_dir

    def _update_status_api(self):
        """
        根据client决定接口api为meta.update_status/meta.tupdate_status
        """
        data = web.input()
        client = data.client if hasattr(data, "client") else web.ctx.env.get('HTTP_CLIENT')
        if client == 'client01':
            return 'meta.update_status'
        else:
            return 'meta.tupdate_status'

# -*- coding: utf-8 -*-
from __future__ import print_function
import web,os
import random
import json
import time
from ..core.basic import Basic
from mainapp.libs.signature import check_sig
from mainapp.models.mongo.metabolome import Metabolome
from mainapp.models.mongo.core.base import Base
from mainapp.models.workflow import Workflow
from biocluster.config import Config


class MetabolomeController(Base):
    def __init__(self, bind_object=None, instant=False):
        super(MetabolomeController, self).__init__(bind_object)
        self._project_type = "metabolome"
        self._instant = instant
        self._post_data = None
        self._sheet_data = None
        self._return_msg = None
        # 下面的Metabolome成了实例对象，拥有多种和数据库交互的函数。
        self.Metabolome = Metabolome(bind_object=bind_object)

    @property
    def data(self):
        return self._post_data

    @property
    def return_msg(self):
        """
        获取Post数据
        """
        return self._return_msg

    @property
    def instant(self):
        """
        任务是否是即时计算
        """
        return self._instant

    @property
    def sheet_data(self):
        """
        获取运行流程所需的Json数据
        :return:
        """
        return self._sheet_data

    def _update_status_api(self):
        """
        根据client决定接口api为
        1. mbio.api.web.metabolome.update_status
        2. mbio.api.web.metabolome.tupdate_status
        """
        data = web.input()
        client = data.client if hasattr(data, "client") else web.ctx.env.get('HTTP_CLIENT')
        if client == 'client01':
            return 'metabolome.update_status'
        else:
            return 'metabolome.tupdate_status'

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
                self.Metabolome.update_status_failed(update_info[i], i)
                print("INFO: 更新主表状态为failed成功: coll:{} _id:{}".format(update_info[i], i))
        except Exception as e:
            print('ERROR:尝试回滚主表状态为failed 失败:{}'.format(e))

    @check_sig
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
        except Exception as e:
            self.roll_back()
            return {"success": False, "info": "运行出错: {}".format(e)}

    def set_sheet_data(self, name, options, main_table_name, task_id, project_sn,
                       module_type="workflow", params=None, to_file=None):
        """
        设置运行所需的Json文档

        :param name: workflow/module/tool相对路径
        :param module_type: workflow/module/tool
        :param main_table_name: 交互分析项主表名称
        :param options: workflow/module/tool参数
        :param params: 交互分析主表params字段
        :param to_file: workflow/module/tool mongo数据转文件
        :return:
        """
        self._post_data = web.input()
        new_task_id = self.get_new_id(task_id)
        self._sheet_data = {
            'id': new_task_id,
            'interaction': True,
            'name': name,  # 需要配置
            'type': module_type,  # 可以配置
            'client': self.data.client,
            'output': self._create_output_dir(task_id, main_table_name),
            'project_sn': project_sn,
            'IMPORT_REPORT_DATA': True,
            'UPDATE_STATUS_API': self._update_status_api(),
            'db_type': 'metabolome',  # 需要根据自己连接的哪个数据库来修改实际的名字
            'options': options  # 需要配置
        }
        if self.instant:
            self._sheet_data["instant"] = True
        if params:
            self._sheet_data["params"] = params
        if to_file:
            self._sheet_data["to_file"] = to_file
        self.workflow_id = new_task_id
        return self._sheet_data

    def _create_output_dir(self, task_id, main_table_name):
        data = web.input()
        task_info = self.Metabolome.get_task_info(task_id)
        client = data.client if hasattr(data, "client") else web.ctx.env.get('HTTP_CLIENT')
        if client == 'client01':
            target_dir = 'sanger'
        else:
            target_dir = 'tsanger'
        '''
        part_dir = str(task_info['member_id']) + \
                   '/' + str(task_info['project_sn']) + '/' + \
                   task_id + '/interaction_results/' + main_table_name
        target_dir = 's3://rerewrweset/files/' + part_dir + '/'
        '''
        part_dir = 'files/' + str(task_info['member_id']) + \
                   '/' + str(task_info['project_sn']) + '/' + \
                   task_id + '/interaction_results/' + main_table_name + "/"
        my_bucket = Config().get_project_region_bucket(project_type="metabolome")
        target_dir = os.path.join(my_bucket, part_dir)
        return target_dir

    def _get_file_path(self, path_postfix, task_id, main_table_name):
        my_bucket = Config().get_project_region_bucket(project_type="metabolome")
        target_dir = self._get_file_dir(path_postfix, task_id, main_table_name, path_type="s3")
        target_path = os.path.join(my_bucket, target_dir)
        return target_path

    def _get_file_abs(self, path_postfix, task_id, main_table_name):
        data = web.input()
        # task_info = self.Metabolome.get_task_info(task_id)
        client = data.client if hasattr(data, "client") else web.ctx.env.get('HTTP_CLIENT')
        if client == 'client01':
            sanger_type = 'sanger'
        else:
            sanger_type = 'tsanger'
        sanger_prefix = Config().get_netdata_config(sanger_type)
        target_dir = sanger_prefix[sanger_type + '_path'] + '/' + self._get_file_dir(path_postfix, task_id,
                                                                                     main_table_name)
        return target_dir

    def _get_file_dir(self, path_postfix, task_id, main_table_name, path_type="rere"):
        task_info = self.Metabolome.get_task_info(task_id)
        if path_type == "rere":
            target_dir = 'rerewrweset/files/' + str(task_info['member_id']) + '/' + str(
                task_info['project_sn']) + '/' + task_id + '/interaction_results/' + main_table_name + path_postfix
        elif path_type == "s3":
            target_dir = 'files/' + str(task_info['member_id']) + '/' + str(
                task_info['project_sn']) + '/' + task_id + '/interaction_results/' + main_table_name + path_postfix + "/"
        return target_dir

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

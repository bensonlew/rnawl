# -*- coding: utf-8 -*-
#__author__ = 'gaohao'

import web
from ..core.basic import Basic
from mainapp.libs.signature import check_sig
from mainapp.models.mongo.bac_comparative import BacComparative
from mainapp.models.mongo.core.base import Base
from mainapp.models.workflow import Workflow
from biocluster.config import Config
import datetime
import os,random
import json

class BacComparativeController(Base):
    """
    设置项目基本信息
    """
    def __init__(self, bind_object=None, instant=False):
        super(BacComparativeController, self).__init__(bind_object)
        self._project_type = "bac_comparative"
        self._instant = instant
        self._post_data = None
        self._sheet_data = None
        self._return_msg = None
        self.bac_comparative = BacComparative(bind_object=bind_object)
        self.config = Config()

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
        根据client决定接口api为metagbin.update_status/metagbin.tupdate_status
        :return:
        """
        data = web.input()
        client = data.client if hasattr(data, "client") else web.ctx.env.get('HTTP_CLIENT')
        if client == 'client01':
            return 'bac_comparative.update_status'
        else:
            return 'bac_comparative.tupdate_status'

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
        except Exception, e:
            self.roll_back()
            return {"success": False, "info": "ERROR: %s" % e}

    def set_sheet_data(self, name, options, main_table_name, task_id, project_sn,module_type="workflow", params=None, to_file=None):
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
            'db_type': 'bac_comparative',  # 需要根据自己连接的哪个数据库来修改实际的名字
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
        """
        创建结果目录
        :param task_id:
        :param main_table_name:
        :return:
        """
        data = web.input()
        task_info = self.bac_comparative.get_task_info(task_id)
        part_dir = 'files/' + str(task_info['member_id']) + \
                      '/' + str(task_info['project_sn']) + '/' + \
                      task_id + '/interaction_results/' + main_table_name + "/"
        my_bucket = Config().get_project_region_bucket(project_type="bac_comparative")
        target_dir = os.path.join(my_bucket, part_dir)
        return target_dir

    def get_new_id(self, task_id):
        """
        按照时间去生成task_id
        :params task_id: 主任务id
        :return:
        """
        new_id = "{}_{}_{}".format(task_id, datetime.datetime.now().strftime("%m%d%H%M%S%f"), random.randint(1, 10000))
        workflow_module = Workflow()
        workflow_data = workflow_module.get_by_workflow_id(new_id)
        if len(workflow_data) > 0:
            return self.get_new_id(task_id)
        print new_id
        return new_id

    def set_sheet_data2(self, name, options, main_table_name, task_id, project_sn, module_type="workflow", params=None,
                       to_file=None):
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
            'name': name,  # 需要配置
            'type': module_type,  # 可以配置
            'client': self.data.client,
            'output': self._create_output_dir2(task_id, main_table_name),
            'project_sn': project_sn,
            'IMPORT_REPORT_DATA': True,
            'UPDATE_STATUS_API': self._update_status_api(),
            'db_type': 'bac_comparative',  # 需要根据自己连接的哪个数据库来修改实际的名字
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

    def _create_output_dir2(self, task_id, main_table_name):
        """
        创建结果目录
        :param task_id:
        :param main_table_name:
        :return:
        """
        data = web.input()
        part_dir = 'files/' + str(data.member_id) + \
                   '/' + str(data.project_sn) + '/' + \
                   task_id + '/interaction_results/' + main_table_name + "/"
        my_bucket = Config().get_project_region_bucket(project_type="bac_comparative")
        target_dir = os.path.join(my_bucket, part_dir)
        return target_dir
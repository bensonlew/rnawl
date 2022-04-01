# -*- coding: utf-8 -*-

import web
from ..core.basic import Basic
from mainapp.libs.signature import check_sig
from mainapp.models.mongo.arghub import Arghub
from mainapp.models.mongo.core.base import Base
from mainapp.models.workflow import Workflow
from biocluster.config import Config
import datetime
import os,random
import json

class ArghubController(Base):
    """
    设置项目基本信息
    """
    def __init__(self, bind_object=None, instant=False):
        super(ArghubController, self).__init__(bind_object)
        self._project_type = "arghub"
        self._instant = instant
        self._post_data = None
        self._sheet_data = None
        self._return_msg = None
        self.arghub = Arghub(bind_object=bind_object)
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
            return 'arghub.update_status'
        else:
            return 'arghub.tupdate_status'

    @check_sig
    def POST(self):
        workflow_client = Basic(data=self.sheet_data, instant=self.instant)
        try:
            run_info = workflow_client.run()
            self._return_msg = workflow_client.return_msg
            return run_info
        except Exception, e:
            self.roll_back()
            return {"success": False, "info": "ERROR: %s" % e}

    def set_sheet_data(self, name, options, main_table_name, task_id, module_type="workflow", params=None, to_file=None, mem_id=None):
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
        self._sheet_data = {
            'id': task_id,
            'name': name,  # 需要配置
            'type': module_type,  # 可以配置
            'client': self.data.client,
            'output': self._create_output_dir(task_id, main_table_name, mem_id),
            'IMPORT_REPORT_DATA': True,
            'UPDATE_STATUS_API': self._update_status_api(),
            'db_type': 'arghub',  # 需要根据自己连接的哪个数据库来修改实际的名字
            'options': options  # 需要配置
        }
        if self.instant:
            self._sheet_data["instant"] = True
        if params:
            self._sheet_data["params"] = params
        if to_file:
            self._sheet_data["to_file"] = to_file
        self.workflow_id = task_id
        return self._sheet_data

    def _create_output_dir(self, task_id, main_table_name, mem_id):
        """
        创建结果目录
        :param task_id:
        :param main_table_name:
        :return:
        """
        data = web.input()
        part_dir = 'files/' + str(mem_id) +\
                      '/' + task_id + '/' +  main_table_name + "/"
        my_bucket = Config().get_project_region_bucket(project_type="arghub")
        target_dir = os.path.join(my_bucket, part_dir)
        return target_dir

    def get_new_id(self, task_id):
        """
        按照时间去生成task_id
        :params task_id: 主任务id
        :return:
        """
        new_id = "{}_{}_{}".format(task_id, datetime.datetime.now().strftime("%m%d%H%M%S"), random.randint(1, 10000))
        workflow_module = Workflow()
        workflow_data = workflow_module.get_by_workflow_id(new_id)
        if len(workflow_data) > 0:
            return self.get_new_id(task_id)
        print new_id
        return new_id


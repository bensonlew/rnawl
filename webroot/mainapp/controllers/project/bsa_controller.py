# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# modified 2018.02.23

import os
import web
import json
import random
import datetime
from ..core.basic import Basic
from biocluster.config import Config
from mainapp.models.workflow import Workflow
from mainapp.libs.signature import check_sig


class BsaController(object):
    def __init__(self, instant=False):
        self._instant = instant
        self._post_data = None
        self._return_msg = None
        self._sheet_data = None

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
    def sheet_data(self):
        """
        获取运行流程所需的Json数据
        :return:
        """
        return self._sheet_data

    @property
    def instant(self):
        """
        任务是否是即时计算
        :return: bool
        """
        return self._instant

    def _update_status_api(self):
        """
        根据client决定接口api为bsa.update_status/bsa.tupdate_status
        :return:
        """
        data = web.input()
        client = data.client if hasattr(data, "client") else web.ctx.env.get("HTTP_CLIENT")
        if client == "client01":
            return "bsa.update_status"
        else:
            return "bsa.tupdate_status"

    @check_sig
    def POST(self):
        workflow_client = Basic(data=self.sheet_data, instant=self.instant)
        try:
            run_info = workflow_client.run()
            self._return_msg = workflow_client.return_msg
            return run_info
        except Exception, e:
            return json.dumps({"success": False, "info": "运行出错：{}".format(e)})

    def set_sheet_data(self, name, member_id, project_sn, task_id, main_table_name, options, module_type="workflow", params=None, to_file=None):
        """
        设定运行参数
        """
        self._post_data = web.input()
        self._sheet_data = {
            "id": self.get_new_id(task_id),
            "name": name,
            "type": module_type,
            'interaction': True,
            "client": self.data.client,
            "IMPORT_REPORT_DATA": True,
            "UPDATE_STATUS_API": self._update_status_api(),
            "db_type": "bsa",
            "options": options,
            "output": self._create_output_dir(member_id, project_sn, task_id, main_table_name, project_type="bsa")
        }
        if self.instant:
            self._sheet_data["instant"] = True
        if params:
            self._sheet_data["params"] = params
        if to_file:
            self._sheet_data["to_file"] = to_file
        return self._sheet_data

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
            return self.get_new_id()
        print new_id
        return new_id

    def _create_output_dir(self, member_id, project_sn, task_id, main_table_name, project_type=None):
        """
        设置接口的上传路径--target_type=all同时上传mnt与s3路径，为s3的时候，直接上传对象存储，为ilustre，只传到我们本地机器
        :param member_id:
        :param project_sn:
        :param task_id:
        :param main_table_name:
        :param target_type:
        :return:
        """
        data = web.input()
        config = Config()
        if config.RGW_ENABLE:
            if project_type:
                bucket = config.get_project_region_bucket(project_type)
            else:
                bucket = config.get_project_region_bucket("main")
            target_dir = os.path.join(bucket, 'files', member_id, project_sn, task_id,
                                      'interaction_results', main_table_name)
        else:
            client = data.client if hasattr(data, "client") else web.ctx.env.get('HTTP_CLIENT')
            if client == 'client01':
                target_dir = 'sanger'
            else:
                target_dir = 'tsanger'
            target_dir += ':rerewrweset/files/' + member_id + '/' + project_sn \
                          + '/' + task_id + '/interaction_results/' + main_table_name
        return target_dir

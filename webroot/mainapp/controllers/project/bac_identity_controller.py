# -*- coding: utf-8 -*-
# __author__ = 'wangzhaoyue'

import web
from ..core.basic import Basic
from mainapp.libs.signature import check_sig
from mainapp.models.workflow import Workflow
import json
import random
import datetime
from biocluster.config import Config


class BacidentityController(object):
    def __init__(self, instant=False):
        self._instant = instant
        self._post_data = None
        self._sheet_data = None
        self._return_msg = None

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

    def _update_status_api(self):
        """
        根据client决定接口api为datasplit.update_status/datasplit.tupdate_status
        """
        data = web.input()
        client = data.client if hasattr(data, "client") else web.ctx.env.get('HTTP_CLIENT')
        if client == 'client01':
            # return 'datasplit.update_status'
            return 'tool_lab.update_status'
        else:
            # return 'datasplit.tupdate_status'
            return 'tool_lab.tupdate_status'

    @check_sig
    def POST(self):
        workflow_client = Basic(data=self.sheet_data, instant=self.instant)
        try:
            run_info = workflow_client.run()
            self._return_msg = workflow_client.return_msg
            return run_info
        except Exception, e:
            return json.dumps({"success": False, "info": "运行出错：{}".format(e)})

    def set_sheet_data(self, name, options, table_id, to_file=None, seq_number=None, module_type="workflow", params=None):
        self._post_data = web.input()
        new_task_id = self.get_new_id(table_id)
        self._sheet_data = {
            'id': new_task_id,
            'stage_id': 0,
            'name': name,
            'type': module_type,
            "USE_DB": True,
            'IMPORT_REPORT_DATA': True,
            'client': self.data.client,
            'IMPORT_REPORT_AFTER_END': False,
            'UPDATE_STATUS_API': self._update_status_api(),
            'instant': False,
            'options': options
        }
        if seq_number:
            self._sheet_data["output"] = self._create_output_dir(seq_number, new_task_id)
        if to_file:
            self._sheet_data["to_file"] = to_file
        if self.instant:
            self._sheet_data['instant'] = True
        if params:
            self._sheet_data['params'] = params
        return self._sheet_data

    def get_new_id(self, task_id):
        new_id = "{}_{}".format(task_id, datetime.datetime.now().strftime("%Y%m%d_%H%M%S"))
        workflow_module = Workflow()
        workflow_data = workflow_module.get_by_workflow_id(new_id)
        if len(workflow_data) > 0:
            return self.get_new_id(task_id)
        return new_id

    def _create_output_dir(self, seq_number, new_task_id):
        """
        设置高通量数据拆分上传文件目录
        """
        # data = web.input()
        # client = data.client if hasattr(data, "client") else web.ctx.env.get('HTTP_CLIENT')
        # if client == 'client01':
        #     target_dir = 's3://datasplit/' + datetime.datetime.now().strftime("%Y") + "/"
        # + seq_number + '/' + new_task_id + '/'
        # else:
        #     target_dir = 's3://rerewrweset/files/datasplit/' + datetime.datetime.now().strftime("%Y") +
        # "/" + seq_number + '/' + new_task_id + '/'
        config = Config()
        bucket = config.get_project_region_bucket(project_type="datasplit").rstrip("/")
        target_dir = '{}/'.format(bucket) + datetime.datetime.now().strftime("%Y") + "/" + seq_number + '/' \
                     + new_task_id + '/'
        return target_dir

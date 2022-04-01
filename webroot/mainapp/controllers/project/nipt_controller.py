# -*- coding: utf-8 -*-
# __author__ = 'hongdong.xuan'

import web
from ..core.basic import Basic
from mainapp.libs.signature import check_sig
from mainapp.models.workflow import Workflow
import datetime
import json
import random


class NiptController(object):
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

    def _update_status_api(self, db):
        """
        根据client决定接口api为denovo.update_status/denovo.tupdate_status
        :param db:
        :return:
        """
        data = web.input()
        client = data.client if hasattr(data, "client") else web.ctx.env.get('HTTP_CLIENT')
        if client == 'client01':
            return '{}.med_report_update'.format(db)
        else:
            return '{}.med_report_tupdate'.format(db)

    @check_sig
    def POST(self):
        workflow_client = Basic(data=self.sheet_data, instant=self.instant)
        try:
            run_info = workflow_client.run()
            self._return_msg = workflow_client.return_msg
            return run_info
        except Exception, e:
            return json.dumps({"success": False, "info": "运行出错：{}".format(e)})

    def set_sheet_data(self, name, options, module_type="workflow", params=None, db_type=None):
        self._post_data = web.input()
        if db_type:
            db = db_type
        else:
            db = 'nipt'
        self._sheet_data = {
            'id': self.get_new_id(),
            'name': name,
            'type': module_type,
            'interaction': True,
            'client': self.data.client,
            'IMPORT_REPORT_DATA': True,
            'UPDATE_STATUS_API': self._update_status_api(db),
            'instant': False,
            'db_type': db,
            'options': options
        }
        print self._sheet_data
        if self.instant:
            self._sheet_data['instant'] = True
        if params:
            self._sheet_data['params'] = params
        return self._sheet_data

    def get_new_id(self):
        """
        按照时间去生成task_id
        :return:
        """
        new_id = 'nipt_{}_{}_{}'.format(datetime.datetime.now().strftime("%m%d%H%M%S"),
                                        random.randint(1000, 10000), random.randint(1, 10000))
        workflow_module = Workflow()
        workflow_data = workflow_module.get_by_workflow_id(new_id)
        if len(workflow_data) > 0:
            return self.get_new_id()
        print new_id
        return new_id

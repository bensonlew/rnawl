# -*- coding: utf-8 -*-
# __author__ = 'hongdongxuan'

from mainapp.models.workflow import Workflow
from ..core.basic import Basic
import datetime
import random
import json
import web


class PtController(object):
    def __init__(self, instant=False):
        """
        lasted modified by hongdong@20180821
        :param instant:
        """
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
        根据client决定接口api为pt_v2.med_report_update/pt_v2.med_report_tupdate
        """
        data = web.input()
        client = data.client if hasattr(data, "client") else web.ctx.env.get('HTTP_CLIENT')
        if client == 'client01':
            return '{}.med_report_update'.format(db)
        else:
            return '{}.med_report_tupdate'.format(db)

    # @check_sig
    def POST(self):
        workflow_client = Basic(data=self.sheet_data, instant=self.instant)
        try:
            run_info = workflow_client.run()
            self._return_msg = workflow_client.return_msg
            return run_info
        except Exception, e:
            return json.dumps({"success": False, "info": "运行出错：{}".format(e)})

    def set_sheet_data_(self, name, options, module_type="workflow", params=None, db_type=None, analysis_name=None):
        """
        设定运行参数
        :param name:
        :param options:
        :param module_type:
        :param params:
        :param db_type:亲子目前的数据库，分为pt与pt_v2，该值其实就是我们的collection要往那个数据库中写
        :param analysis_name:亲子分析的接口可以单独进行call snp或者进行家系分析，只针对V2版，该信息会用到结果目录中，用于区分分析
        :return:
        """
        self._post_data = web.input()
        if db_type:
            db = db_type
        else:
            db = 'pt'
        self._sheet_data = {
            'id': self.get_new_id(analysis_name) if analysis_name else self.get_new_id(),
            'name': name,
            'type': module_type,
            'client': self.data.client,
            'interaction': True,
            'IMPORT_REPORT_DATA': True,
            'UPDATE_STATUS_API': self._update_status_api(db),
            'instant': False,
            'db_type': db,
            'options': options,
            'output': self._create_output_dir(name)
        }
        if self.instant:
            self._sheet_data['instant'] = True
        if params:
            self._sheet_data['params'] = params
        return self._sheet_data

    def get_new_id(self, analysis_name=None):
        """
        按照时间去生成task_id
        :param analysis_name:
        :return:
        """
        if analysis_name:
            new_id = '{}_{}_{}_{}'.format(analysis_name, datetime.datetime.now().strftime("%m%d%H%M%S"),
                                          random.randint(1000, 10000), random.randint(1, 10000))
        else:
            new_id = 'pt_{}_{}_{}'.format(datetime.datetime.now().strftime("%m%d%H%M%S"),
                                          random.randint(1000, 10000), random.randint(1, 10000))
        workflow_module = Workflow()
        workflow_data = workflow_module.get_by_workflow_id(new_id)
        if len(workflow_data) > 0:
            return self.get_new_id()
        print new_id
        return new_id

    def _create_output_dir(self, name):
        """
        设置医学上传文件目录
        add by yuguo

        """
        data = web.input()
        client = data.client if hasattr(data, "client") else web.ctx.env.get('HTTP_CLIENT')
        if client == 'client01':
            target_dir = 'sanger'
        else:
            target_dir = 'tsanger'
        target_dir += ':rerewrweset/MEDfiles/' + name
        return target_dir

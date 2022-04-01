# -*- coding: utf-8 -*-
# __author__ = 'hongdongxuan'
# last modified by hd @ 20180408

from mainapp.models.workflow import Workflow
from biocluster.config import Config
from ..core.basic import Basic
import datetime
import random
import json
import web
import os


class DnaController(object):
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
        根据client决定接口api为update_status/tupdate_status
        """
        data = web.input()
        client = data.client if hasattr(data, "client") else web.ctx.env.get('HTTP_CLIENT')
        if client == 'client01':
            return '{}.update_status'.format(db)
        else:
            return '{}.tupdate_status'.format(db)

    # @check_sig
    def POST(self):
        workflow_client = Basic(data=self.sheet_data, instant=self.instant)
        try:
            run_info = workflow_client.run()
            self._return_msg = workflow_client.return_msg
            return run_info
        except Exception, e:
            return json.dumps({"success": False, "info": "运行出错：{}".format(e)})

    def set_sheet_data(self, name, member_id, project_sn, task_id, main_table_name, options, module_type="workflow",
                       params=None, db_type=None, analysis_name=None, to_file=None, target_output=False):
        """
        设定运行参数
        :param name:
        :param options:
        :param module_type:
        :param params:
        :param member_id:
        :param project_sn:
        :param task_id:
        :param main_table_name:
        :param to_file:
        :param db_type:目前只有wgs，要注意mbio.web下面的项目文件夹名字 要去mongo的连接一致
        :param analysis_name:区分是wgs还是遗传图谱
        :param target_output:当接口直接调用tool或者module的时候，结果文件要传到指定的路径下，但是tool or module中获取不到
        该路径，当该参数为True的时候，会在options添加tool中文件要上传的路径，默认是不需要的。
        :return:
        """
        self._post_data = web.input()
        if db_type:
            db = db_type
        else:
            db = 'dna_wgs'
        self._sheet_data = {
            'id': self.get_new_id(task_id, analysis_name) if analysis_name else self.get_new_id(task_id),
            'name': name,
            'interaction': True,
            'type': module_type,
            'client': self.data.client,
            'IMPORT_REPORT_DATA': True,
            'UPDATE_STATUS_API': self._update_status_api(db),
            'instant': False,
            'db_type': db,
            'options': options,
            'output': self._create_output_dir(member_id, project_sn, task_id, main_table_name, project_type=db)
        }
        if self.instant:
            self._sheet_data['instant'] = True
        if params:
            self._sheet_data['params'] = params
        if to_file:
            self._sheet_data["to_file"] = to_file
        if target_output:
            self._sheet_data['options'].update({"target_path": self._sheet_data['output']})
        return self._sheet_data

    def get_new_id(self, task_id, analysis_name=None):
        """
        按照时间去生成task_id
        :param analysis_name:
        :return:
        """
        if analysis_name:
            new_id = '{}_{}_{}'.format(task_id, analysis_name+datetime.datetime.now().strftime("%m%d%H%M%S%f"),
                                       random.randint(1000, 10000))
        else:
            new_id = '{}_{}_{}'.format(task_id, datetime.datetime.now().strftime("%m%d%H%M%S%f"),
                                       random.randint(1000, 10000))
        workflow_module = Workflow()
        workflow_data = workflow_module.get_by_workflow_id(new_id)
        if len(workflow_data) > 0:
            if analysis_name:
                return self.get_new_id(task_id, analysis_name)
            else:
                return self.get_new_id(task_id)
        print new_id
        return new_id

    def _create_output_dir(self, member_id, project_sn, task_id, main_table_name, project_type=None):
        """
        设置接口的上传路径--target_type=all同时上传mnt与s3路径，为s3的时候，直接上传对象存储，为ilustre，只传到我们本地机器
        :param member_id:
        :param project_sn:
        :param task_id:
        :param main_table_name:
        :param project_type:
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
            target_dir += ':rerewrweset/files/' + member_id + '/' + project_sn\
                          + '/' + task_id + '/interaction_results/' + main_table_name
        return target_dir

# -*- coding: utf-8 -*-
# __author__ = 'HD'
# last modified by hd @ 20200402
# last modified by hd @ 20200520 增加即时分析，当即返回文件路径信息

from mainapp.models.workflow import Workflow
# from biocluster.config import Config
from mainapp.config.db import Config
# from ..core.basic import Basic
from mainapp.controllers.core.basic import Basic
import datetime
import random
import json
import web
import os
import time
import re


class ToollabController(object):
    def __init__(self):
        self.instant = True
        self._post_data = None
        self._sheet_data = None
        self._return_msg = None
        self.client = 'client01'
        self.id = ''
        self.name = ''
        self.type = 'workflow'
        self.project_sn = ''

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

    def _update_status_api(self, db):
        """
        根据client决定接口api为update_status/tupdate_status
        """
        if self.client == 'client01':
            return '{}.update_status'.format(db)
        else:
            return '{}.tupdate_status'.format(db)

    # @check_sig
    def POST(self):
        workflow_client = Basic(data=self.sheet_data, instant=self.instant)
        try:
            run_info = workflow_client.run()
            if run_info['success'] and self.instant:
                post_data = self.get_post_data()
                run_info['post_data'] = post_data
            self._return_msg = workflow_client.return_msg
            return run_info
        except Exception, e:
            return {"success": False, "info": "运行出错：{}".format(e)}

    def get_post_data(self):
        """
        小工具即时分析的时候，需要将运行完成的结果文件信息，返还给前端进行插入mysql中。
        :return:
        """
        work_dir = self.worker_dir()
        timestr = str(time.strftime('%Y%m%d', time.localtime(time.time())))
        if self.type == 'workflow':
            # 大驼峰命令
            module_name = "".join(' '.join(self.name.strip().split('.')[-1].split('_')).title().split(' '))
        else:
            module_name = 'Single'
        work_dir = work_dir + "/" + timestr + "/" + module_name + "_" + self.id
        print work_dir
        if os.path.exists(work_dir):
            post_result_data = os.path.join(work_dir, 'post_result_data.json')
            if os.path.exists(post_result_data):
                with open(post_result_data, 'r') as r:
                    data = r.read()
            else:
                return self.make_no_file_post_data()
        else:
            return
        return self.make_new_post_data(json.loads(data))

    def make_new_post_data(self, data):
        """
        按照前端要求组建正确格式的数据结构
        :return:
        """
        data = data[u'data']
        post_data = dict()
        content = {
            "basis": {
                "task_id": self._get_main_task_id()
            }
        }
        if u"task" in data[u'sync_task_log'].keys():
            format_task = data[u"sync_task_log"][u"task"]
            format_task[u'task_id'] = self._get_main_task_id()
            content["basis"] = format_task
        if u'files' in data[u'sync_task_log'].keys():
            content['files'] = data[u"sync_task_log"][u"files"]
        if u'dirs' in data[u'sync_task_log'].keys():
            content['dirs'] = data[u'sync_task_log'][u'dirs']
        if u'base_path' in data[u'sync_task_log'].keys():
            content['basis']['base_path'] = data[u'sync_task_log'][u'base_path']
        if u"region" in data[u'sync_task_log'].keys():
            content['basis']['region'] = data[u'sync_task_log'][u'region']
        if u'log' in data[u'sync_task_log'].keys():
            content['log'] = data[u'sync_task_log'][u'log']
        content['basis']['project_sn'] = self.project_sn
        post_data['sync_task_log'] = content
        return post_data

    def make_no_file_post_data(self):
        """
        有些小工具没有文件路径，所以不会生成post_result_data.json文件，但是运行成功了，所以这里要返回状态信息
        {\"sync_task_log\": {\"log\": [{\"status\": \"finish\", \"time\": \"2020-06-16 17:37:34\",
        \"desc\": \"Job has been finished\"}], \"basis\": {\"status\": \"finish\",
        \"created_ts\": \"2020-06-16 17:37:34\", \"task_id\": \"47a4qshpi7ph7o7d4cco78cibb\"}}}
        :return:
        """
        post_data = dict()
        content = {
            "basis": {
                "task_id": self._get_main_task_id(),
                'status': 'finish',
                'created_ts': datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            },
            'log': [
                {
                    'status': 'finish',
                    'time': datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                    'desc': 'Job has been finished'
                }
            ]
        }
        post_data['sync_task_log'] = content
        return post_data

    def _get_main_task_id(self):
        my_id = re.split('_', self.id)
        my_id.pop(-1)
        my_id.pop(-1)
        return "_".join(my_id)

    def worker_dir(self):
        # config = Config()
        # return config.WORK_DIR
        return Config().get_work_dir()

    def set_sheet_data(self, name, project_sn, task_id, main_table_name, options, module_type="workflow", params=None,
                       db_type=None, analysis_name=None, to_file=None, client="client01", target_output=False,
                       instant=True):
        """
        设定运行参数
        :param name:
        :param task_id:
        :param options:
        :param main_table_name:
        :param to_file:
        :param module_type:
        :param instant:
        :param params:
        :param client:
        :param project_sn:
        :param db_type:目前默认为toollab，要注意mbio.web下面的项目文件夹名字 要与mongo的连接key一致
        :param analysis_name:区分是那个小工具
        :param target_output:当接口直接调用tool或者module的时候，结果文件要传到指定的路径下，但是tool or module中获取不到
        该路径，当该参数为True的时候，会在options添加tool中文件要上传的路径，默认是不需要的。
        :return:
        """
        self._post_data = web.input()
        self.instant = instant
        self.client = client
        self.name = name
        self.project_sn = project_sn
        if db_type:
            db = db_type
        else:
            db = 'tool_lab'
        self._sheet_data = {
            'id': self.get_new_id(task_id, analysis_name) if analysis_name else self.get_new_id(task_id),
            'name': name,
            'type': module_type,
            'interaction': True,
            'client': self.client,
            'IMPORT_REPORT_DATA': True,
            'UPDATE_STATUS_API': self._update_status_api(db),
            'instant': False,
            'db_type': db,
            "project_sn": project_sn,
            'options': options,
            'output': self._create_output_dir(task_id, project_sn, main_table_name, project_type=db)
        }
        if self.instant:
            self._sheet_data['instant'] = True
        if params:
            self._sheet_data['params'] = params
        if to_file:
            self._sheet_data["to_file"] = to_file
        if target_output:
            self._sheet_data['options'].update({"target_path": self._sheet_data['output']})
        self.id = self._sheet_data['id']
        self.type = self._sheet_data['type']
        # print self._sheet_data
        return self._sheet_data

    def get_new_id(self, task_id, analysis_name=None):
        """
        按照时间去生成task_id
        :param task_id:
        :param analysis_name:
        :return:
        """
        if analysis_name:
            new_id = '{}_{}_{}'.format(task_id, analysis_name, datetime.datetime.now().strftime("%m%d%H%M%S%f"),
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

    def _create_output_dir(self, task_id, project_sn, main_table_name, project_type=None):
        """
        设置接口的上传路径--为s3的时候，直接上传对象存储，为ilustre，只传到我们本地机器
        :param task_id:
        :param main_table_name:
        :param project_type:
        :return:
        """
        config = Config()
        if config.RGW_ENABLE:
            if project_type:
                bucket = config.get_project_region_bucket(project_type)
            else:
                bucket = config.get_project_region_bucket("main")
            target_dir = os.path.join(bucket, 'files', project_sn, task_id)
        else:
            if self.client == 'client01':
                target_dir = 'sanger'
            else:
                target_dir = 'tsanger'
            target_dir += ':rerewrweset/files/' + project_sn + '/' + task_id + '/' + main_table_name
        return target_dir

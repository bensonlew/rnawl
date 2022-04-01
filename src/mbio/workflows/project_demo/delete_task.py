# -*- coding: utf-8 -*-
#__author__ = 'qingchen.zhang'
import shutil
import os
import gevent
from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError
from mbio.packages.project_demo.clean_mongo import CleanMongo

class DeleteTaskWorkflow(Workflow):
    """
    1.用于页面删除交互分析的记录同时删除mongo数据；
    2.用于页面删除task或者批量删除，同时删除mongo数据；
    注意目前只能删除task_id的任务，删除交互分析的功能还未能实现
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(DeleteTaskWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "task_id", "type": 'string', "default": ''},#删除任务的id
            {"name": "project_type", "type": 'string', "default": ''},#删除项目类型
            {"name": "location", "type": 'string', "default":'workflow'},#工作流还是交互分析
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())

    def check(self):
        if not self.option('task_id'):
            raise OptionError("Must set task id of what you want to delate")
        if not self.option('project_type'):
            raise OptionError("Must set project type of what you want to delate")

    def delete_task(self):
        """
        删除整个工作流任务的mongo数据
        :return:
        """
        if os.path.exists(self.work_dir + "/task"):
            shutil.rmtree(self.work_dir + "/task")
        os.mkdir(self.work_dir + "/task")
        project_type = self.option('project_type')
        cm = CleanMongo(None, project_type)
        task_list = self.option('task_id').split(',')
        for task_id in task_list:
            task_report_file = os.path.join(self.work_dir, "task/%s.txt"%("record"))
            cm.get_detail_table_info_from_mongo('table_relation', main_id_name='main_id', main_id_loc='asso_id')
            cm.rm_task(task_id, record_file=task_report_file)
            cm.rm_task(task_id)
        gevent.spawn_later(5, self.end)

    def delete_interactive_task(self):
        """
        删除某个任务下的交互分析的数据
        :return:
        """
        if os.path.exists(self.work_dir + "task"):
            shutil.rmtree(self.work_dir + "task")
        os.mkdir(self.work_dir + "task")
        project_type = self.option('project_type')
        cm = CleanMongo(None, project_type)
        task_list = self.option('task_id').split(',')
        for task_id in task_list:
            main_report_file = os.path.join(self.work_dir, "task/%s_main.txt"%(task_id))
            detail_report_file = os.path.join(self.work_dir, "task/%s_detail.txt"%(task_id))
            cm.run(task_id, conf_db='table_relation', detail_log=detail_report_file, main_log=main_report_file)
        gevent.spawn_later(5, self.end)

    def run(self):
        if self.option('location') == 'workflow':
            self.delete_task()
        else:
            self.delete_interactive_task()
        super(DeleteTaskWorkflow, self).run()


    def end(self):
        super(DeleteTaskWorkflow, self).end()
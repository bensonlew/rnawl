# -*- coding: utf-8 -*-
#__author__ = 'qingchen.zhang'

""" 从文件中获取数据统计信息"""
import os
import datetime
from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError


class DataStatWorkflow(Workflow):
    """
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(DataStatWorkflow,self).__init__(wsheet_object)
        options = [
            {"name": "data_table", "type": "infile", "format": "sequence.profile_table"}
        ]
        #self.workflow_output_tmp = self.sheet.output
       # if re.match(r'tsanger:',self.workflow_output_tmp):
            #self.workflow_output = self.workflow_output_tmp.replace('tsanger:','/mnt/ilustre/tsanger-data')
        #else:
            #self.workflow_output =self.workflow_output_tmp.replace('sanger:','/mnt/ilustre/data')
        #self.project_sn = self.sheet.project_sn
        #self.task_id =self._sheet.id
        self.add_option(options)
        self.set_options(self._sheet.options())
        

    def check_options(self):

        pass

    def run(self):
        #if self._sheet.type == "tool":
        self.stat = self.add_tool("test.data_stat")
        #self._task.sheet = self._sheet
        #self.output_dir = self.task.output_dir
        self.output_dir = self.stat.output_dir
        options = {
            'data_table': self.option('data_table')
        }
        #self._task.set_options(options)
        self.stat.set_options(options)
        self.stat.on('end',self.end)
        self.stat.run()
        super(DataStatWorkflow, self).run()

    def set_db(self):
        """
        将运行结果保存到mongo数据库中
        :return:
        """
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        print self.get_upload_files()
        super(DataStatWorkflow, self).end()



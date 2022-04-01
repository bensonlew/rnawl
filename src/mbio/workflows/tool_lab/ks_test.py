# -*- coding: utf-8 -*-
# __author__ = 'zzg'

import os
from biocluster.workflow import Workflow
import datetime
import unittest
import types
from bson.objectid import ObjectId


class KsTestWorkflow(Workflow):
    """
    正态性分布检验

    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(KsTestWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "input_data", "type": "infile","format": "tool_lab.simple"},
            {"name": "main_id", "type": "string"},
            {'name': 'update_info', 'type': 'string'}
        ]
        self.add_option(options)
        self.revise_infiles()
        self.set_options(self._sheet.options())
        self.ks = self.add_tool("tool_lab.normal_test")

    def check_options(self):    
        return True
    
    def run(self):
        self.run_kstest()
        super(KsTestWorkflow, self).run()
        
    def run_kstest(self):
        opts = {
            "input_data": self.option("input_data")
        }
        self.ks.set_options(opts)
        self.ks.on("end", self.set_output)
        self.ks.run()

    def set_output(self):
        for file in os.listdir(self.ks.output_dir):
            os.link(os.path.join(self.ks.output_dir, file), os.path.join(self.output_dir, file))
        self.set_db()
        
    def set_db(self):
        self.logger.info("开始导表")
        api_ks = self.api.api("tool_lab.ks_test")
        api_ks.add_KsTest(self.output_dir,main_id = self.option('main_id'))
        self.logger.info("导表结束")
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        super(KsTestWorkflow, self).end()
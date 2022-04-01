# -*- coding: utf-8 -*-
# __author__ = 'shicaiping'

import os
import re
from biocluster.workflow import Workflow
import datetime
import unittest
import types
from bson.objectid import ObjectId
from biocluster.core.exceptions import OptionError
from biocluster.config import Config


class SraExplorerWorkflow(Workflow):
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(SraExplorerWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "accession", "type": "string", "default": ""},  # 多个以逗号分隔
            {"name": "keyword", "type": "string", "default": ""},
            {"name": "start", "type": "int", "default": 0},
            {"name": "stop", "type": "int", "default": 20},
            {'name': 'update_info', 'type': 'string'},
            {'name': 'main_id', 'type': 'string'}
        ]
        self.add_option(options)
        self.revise_infiles()
        self.tool = self.add_tool("tool_lab.sra_explorer")
        self.set_options(self._sheet.options())

    def run(self):
        self.run_tool()
        super(SraExplorerWorkflow, self).run()

    def check_options(self):
        if not (self.option("accession") or self.option("keyword")):
            raise OptionError("必填参数输入为空")
        return True

    def run_tool(self):
        options = {
            "accession": self.option("accession"),
            "keyword": self.option("keyword"),
            "start": self.option("start"),
            "stop": self.option("stop"),
        }
        self.tool.set_options(options)
        self.tool.on('end', self.set_output)
        self.tool.run()

    def set_output(self):
        for file in os.listdir(self.tool.output_dir):
            os.link(os.path.join(self.tool.output_dir, file), os.path.join(self.output_dir, file))
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "SRA Explorer结果目录", 0],
            ["./*metadata*", "", "SRA Explorer数据文件表", 0,],
        ])
        super(SraExplorerWorkflow, self).end()


class TestFunction(unittest.TestCase):
    '''
    This is test for the workflow. Just run this script to do test.
    '''

    def test(self):
        from mbio.workflows.tool_lab.sra_explorer import SraExplorerWorkflow
        from biocluster.wsheet import Sheet
        import random
        data = {
            "id": "SraExplorer_" + str(random.randint(1, 10000)),
            "type": "workflow",
            "name": "tool_lab.sra_explorer",
            "options": dict(
                accession="SRP008280,SRP008331,SRP010280",
            )
        }
        wsheet = Sheet(data=data)
        wf =SraExplorerWorkflow(wsheet)
        wf.sheet.id = 'sra_explorer'
        wf.sheet.project_sn = 'sra_explorer'
        wf.IMPORT_REPORT_DATA = False
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)

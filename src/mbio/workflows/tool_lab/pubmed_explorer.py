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


class PubmedExplorerWorkflow(Workflow):
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(PubmedExplorerWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "gene_str", "type": "string", "default": ""},  # 多个以逗号分隔
            {"name": "gene_file", "type": "infile", "format": "small_rna.common"},
            {"name": "type", "type": "string", "default": ""},  # str or file
            {"name": "date1", "type": "string", "default": ""},
            {"name": "date2", "type": "string", "default": ""},
            {"name": "if1", "type": "string", "default": ""},
            {"name": "if2", "type": "string", "default": ""},
            {'name': 'update_info', 'type': 'string'},
            {'name': 'main_id', 'type': 'string'}
        ]
        self.add_option(options)
        self.revise_infiles()
        self.tool = self.add_tool("tool_lab.pubmed_explorer")
        self.set_options(self._sheet.options())

    def run(self):
        self.run_tool()
        super(PubmedExplorerWorkflow, self).run()

    def check_options(self):
        if self.option("type") == "str":
            if not self.option("gene_str"):
                raise OptionError("请手动输入基因名称")
        else:
            if not self.option("gene_file").is_set:
                raise OptionError("请上传基因名称列表")
        return True

    def run_tool(self):
        options = {
            "type": self.option("type"),
            "date1": self.option("date1"),
            "date2": self.option("date2"),
            "if1": self.option("if1"),
            "if2": self.option("if2"),
        }
        if self.option("type") == "str":
            if "\n" in self.option("gene_str"):
                tmp_list = ",".join(self.option("gene_str").strip().split("\n"))
                for i in tmp_list:
                    if i == "":
                        tmp_list.remove(i)
                gene_str = str(tmp_list)
                options.update({"gene_str": gene_str})
            else:
                options.update({"gene_str": self.option("gene_str")})
        else:
            options.update({"gene_file": self.option("gene_file")})
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
            ["./pubmed.xlsx", "", "基因文献查询结果表", 0],
        ])
        super(PubmedExplorerWorkflow, self).end()


class TestFunction(unittest.TestCase):
    '''
    This is test for the workflow. Just run this script to do test.
    '''

    def test(self):
        from mbio.workflows.tool_lab.pubmed_explorer import PubmedExplorerWorkflow
        from biocluster.wsheet import Sheet
        import random
        data = {
            "id": "PubmedExplorer_" + str(random.randint(1, 10000)),
            "type": "workflow",
            "name": "tool_lab.pubmed_explorer",
            "options": dict(
                gene_file="/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/pubmed_explorer/gene.txt",
                type="file",
                date1="2020",
                date2="2021",
                if1="3",
                if2="10"
            )
        }
        wsheet = Sheet(data=data)
        wf =PubmedExplorerWorkflow(wsheet)
        wf.sheet.id = 'pubmed_explorer'
        wf.sheet.project_sn = 'pubmed_explorer'
        wf.IMPORT_REPORT_DATA = False
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)

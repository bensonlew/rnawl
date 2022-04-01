# -*- coding: utf-8 -*-
# __author__ = 'zoujiaxun'

"""表格筛选"""

import os
from biocluster.workflow import Workflow
import datetime
import unittest
import types
from bson.objectid import ObjectId
import pandas as pd

class HclusterWorkflow(Workflow):
    """
    处理表格
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(HclusterWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "method", "type": "string", "default": "bray_curtis"},
            {"name": "otutable", "type": "infile",
             "format": "meta.otu.otu_table"},
            {"name": "linkage", "type": 'string', "default": "average"},
            {"name": "group", "type": 'infile', "format": "ref_rna_v2.common"},
            {"name": "sep", "type": 'string', "default": "tab"},
            {'name': 'main_title', 'type': 'string'},
            {'name': 'main_id', 'type': 'string'},
            {'name': 'update_info', 'type': 'string'}
        ]
        self.add_option(options)
        self.revise_infiles()
        # self.tool = self.add_tool("tool_lab.table_standard")
        self.set_options(self._sheet.options())

    def run(self):
        self.run_hcluster()
        super(HclusterWorkflow, self).run()

    def run_hcluster(self):

        self.hcluster = self.add_module('tool_lab.hcluster')

        self.hcluster.set_options({
            "method": self.option("method").lower(),
            "otutable": self.option("otutable"),
            "linkage": self.option('linkage').lower(),
            "sep": self.option('sep')
        })
        self.hcluster.on('end', self.set_db)
        self.hcluster.run()

    def set_db(self):
        """
        保存结果表到mongo数据库中
        """

        hcluster = self.api.api("tool_lab.hcluster")
        # add result info
        hcluster_tre = os.path.join(self.hcluster.output_dir, 'hcluster.tre')
        hcluster.add_hcluster(hcluster_tre, self.option('group').path, main_id=self.option('main_id'), main_title=self.option('main_title'))
        self.set_output()
        self.end()
    # def set_db(self):
    #     """
    #     保存结果标准化数据到mongo数据库中
    #     """
    #     # record_id = self.option("standard_id")
    #     # if isinstance(record_id, types.StringTypes):
    #     #     record_id = ObjectId(record_id)
    #     # elif isinstance(record_id, ObjectId):
    #     #     record_id = record_id
    #     # else:
    #     #     self.set_error("main_id参数必须为字符串或者ObjectId类型!", code="13701601")
    #     standard = self.api.api("tool_lab.table_kit")
    #     # add result info
    #     print self.option('transposed_table').path
    #     standard.add_table_standard()
    #     self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
             [".", "", "层级聚类文件",0],
             ["hcluster.tre", "tre", "层级聚类树文件",0],
             ["otutable.txt", "xls", "样本距离矩阵文件", 0]
        ])
        super(HclusterWorkflow, self).end()

    def set_output(self):
        for file in os.listdir(self.hcluster.output_dir):
            os.link(os.path.join(self.hcluster.output_dir, file), os.path.join(self.output_dir, file))
class TestFunction(unittest.TestCase):
    '''
    This is test for the workflow. Just run this script to do test.
    '''

    def test(self):
        from mbio.workflows.tool_lab.hcluster import HclusterWorkflow
        from biocluster.wsheet import Sheet
        import random
        data = {
            'id': 'hcluster_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'workflow',
            'name': 'tool_lab.hcluster',
            'options': {
                "method" : "bray_curtis",
                "otutable" : "/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/otutable.txt",
                "linkage" : "average",
                'group':"/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/GROUP.txt",
                "sep": 'tab'
            }
        }
        wsheet = Sheet(data=data)
        wf =HclusterWorkflow(wsheet)
        wf.sheet.id = 'cluster'
        wf.sheet.project_sn = 'cluster'
        wf.IMPORT_REPORT_DATA = False
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)

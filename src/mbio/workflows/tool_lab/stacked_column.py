# -*- coding: utf-8 -*-
# __author__ = 'fwy'

import os
from biocluster.workflow import Workflow
import datetime
import unittest
import types
from bson.objectid import ObjectId


class StackedColumnWorkflow(Workflow):
    """
    Used for fasta file seq stat .

    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(StackedColumnWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "raw_table", "type": "infile", "format": "ref_rna_v2.common"},  # fasta文件
            {"name": "combine_value", "type": "float", "default": 0.05},  # 过滤阈值，低于该阈值的归为others
            {"name": "sep", "type": "string", "default": "tab"},  # 表格分隔符
            {"name": "main_id" ,"type": "string"},
            # {"name": "min_len", "type": "int", },
            {'name': 'update_info', 'type': 'string'}
        ]
        self.add_option(options)
        self.revise_infiles()
        self.tool = self.add_tool("tool_lab.stacked_column")
        self.set_options(self._sheet.options())

    def run(self):
        self.run_tool()
        super(StackedColumnWorkflow, self).run()

    def run_tool(self):
        opts = {
            'input_file': self.option("raw_table").prop["path"],
            'filter_value': self.option('combine_value'),
            'sep': self.option('sep')
        }
        self.tool.set_options(opts)
        self.tool.on('end', self.set_db)
        self.tool.run()

    def set_db(self):
        """
                保存结果标准化数据到mongo数据库中
                """
        record_id = self.option("main_id")
        if isinstance(record_id, types.StringTypes):
            record_id = ObjectId(record_id)
        elif isinstance(record_id, ObjectId):
            record_id = record_id
        else:
            self.set_error("main_id参数必须为字符串或者ObjectId类型!", code="13701601")
        stacked = self.api.api("tool_lab.stacked_column")
        # add result info
        stacked_result=os.path.join(self.tool.output_dir,"stacked_column.xls")
        stacked.add_stacked_detail(stacked_result,new_otu_id=record_id)
        self.set_output()


    def set_output(self):
        for file in os.listdir(self.tool.output_dir):
            os.link(os.path.join(self.tool.output_dir, file), os.path.join(self.output_dir, file))
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "堆叠柱形图结果文件",0],
            # [r'.*\.cpm2tpm\.xls', 'xls', '定量指标cpm转tpm结果文件', 0],
            # [r'.*\.count2tpm\.xls', 'xls', '定量指标count转tpm结果文件', 0],
            # [r'.*\.count2cpm\.xls', 'xls', '定量指标count转cpm结果文件', 0],
            # [r'.*\.count2fpkm\.xls', 'xls', '定量指标count转fpkm结果文件', 0],
            # [r'.*\.fpkm2tpm\.xls', 'xls', '定量指标fpkm转tpm结果文件', 0],
            # [r'.*\.cpm2fpkm\.xls', 'xls', '定量指标cpm转fpkm结果文件', 0],
            # [r'.*\.fpkm2cpm\.xls', 'xls', '定量指标fpkm转cpm结果文件', 0],
        ])
        super(StackedColumnWorkflow, self).end()


class TestFunction(unittest.TestCase):
    '''
    This is test for the workflow. Just run this script to do test.
    '''

    def test(self):
        from mbio.workflows.tool_lab.stacked_column import StackedColumnWorkflow
        from biocluster.wsheet import Sheet
        import random
        test_dir = "/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/fungi/Saccharomyces_cerevisiae/Ensemble_release_39/"
        data = {
            "id": "stacked_column" + str(random.randint(1, 10000)),
            "type": "workflow",
            "name": "tool_lab.stacked_column",
            "options": dict(
                input_file="/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/tool_lab/stacked_column/in_otu_table.xls",
                filter_value=0.3,
                sep="tab"
            )
        }
        wsheet = Sheet(data=data)
        wf =StackedColumnWorkflow(wsheet)
        wf.sheet.id = 'stacked'
        wf.sheet.project_sn = 'stacked'
        wf.IMPORT_REPORT_DATA = False
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)

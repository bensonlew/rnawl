# -*- coding: utf-8 -*-
# __author__ = 'fwy'

import os
from biocluster.workflow import Workflow
import datetime
import unittest
import types
from bson.objectid import ObjectId


class FormatConversionWorkflow(Workflow):
    """
    Used for cds to protein code

    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(FormatConversionWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "input_file", "type": "infile", "format": "ref_rna_v2.common"},  # 输入文件，可以是gtf/gff文件
            {"name": "input_format", "type": "string", "default": "gff"},  # ["gtf","gff"]
            {"name": "output_format", "type": "string", "default": "gtf"},  # ["gtf","gff","bed"]
            {"name": "main_id", "type": "string",},
            {'name': 'update_info', 'type': 'string'}
        ]
        self.add_option(options)
        self.revise_infiles()
        self.tool = self.add_tool("tool_lab.format_conversion")
        self.set_options(self._sheet.options())

    def run(self):
        self.run_tool()
        super(FormatConversionWorkflow, self).run()

    def run_tool(self):
        opts = {
            'input_file': self.option('input_file'),
            'input_format': self.option('input_format'),
            'output_format' : self.option("output_format"),

        }
        self.tool.set_options(opts)
        self.tool.on('end', self.set_output)
        self.tool.run()

    def set_output(self):
        for file in os.listdir(self.tool.output_dir):
            os.link(os.path.join(self.tool.output_dir, file), os.path.join(self.output_dir, file))
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        # result_dir.add_relpath_rules([
        #     [".", "", "基因表达量定量指标转换结果文件",0],
        #     [r'.*\.cpm2tpm\.xls', 'xls', '定量指标cpm转tpm结果文件', 0],
        #     [r'.*\.count2tpm\.xls', 'xls', '定量指标count转tpm结果文件', 0],
        #     [r'.*\.count2cpm\.xls', 'xls', '定量指标count转cpm结果文件', 0],
        #     [r'.*\.count2fpkm\.xls', 'xls', '定量指标count转fpkm结果文件', 0],
        #     [r'.*\.fpkm2tpm\.xls', 'xls', '定量指标fpkm转tpm结果文件', 0],
        #     [r'.*\.cpm2fpkm\.xls', 'xls', '定量指标cpm转fpkm结果文件', 0],
        #     [r'.*\.fpkm2cpm\.xls', 'xls', '定量指标fpkm转cpm结果文件', 0],
        # ])
        super(FormatConversionWorkflow, self).end()


class TestFunction(unittest.TestCase):
    '''
    This is test for the workflow. Just run this script to do test.
    '''

    def test(self):
        from mbio.workflows.tool_lab.format_conversion import FormatConversionWorkflow
        from biocluster.wsheet import Sheet
        import random
        test_dir = "/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/fungi/Saccharomyces_cerevisiae/Ensemble_release_39/"
        data = {
            "id": "test" + str(random.randint(1, 10000)),
            "type": "workflow",
            "name": "tool_lab.format_conversion",
            "options": dict(

                output_format="bed",
                input_type="gtf",
                input_file="/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/tool_lab_batch2/gtf_bed_gff/data/gtf/testsmall.gtf",
            )
        }
        wsheet = Sheet(data=data)
        wf =FormatConversionWorkflow(wsheet)
        wf.sheet.id = 'format_conversion'
        wf.sheet.project_sn = 'format_conversion'
        wf.IMPORT_REPORT_DATA = False
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)

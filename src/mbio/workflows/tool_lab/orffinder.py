# -*- coding: utf-8 -*-
# __author__ = 'shicaiping'

import os
from biocluster.workflow import Workflow
import datetime
import unittest
import types
from bson.objectid import ObjectId
from biocluster.core.exceptions import OptionError


class OrffinderWorkflow(Workflow):
    """
    在目标DNA序列中搜寻开发阅读框，可以输出每个ORF所在的区域，并翻译成对应的蛋白序列
    此工具可以为新预测的DNA序列查找潜在的蛋白编码区
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(OrffinderWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "input_file", "type": "infile", "format": "sequence.fasta"},  # FASTA序列文件
            {"name": "start_codon", "type": "int", "default": 2},
            # 0 = "ATG" only
            # 1 = "ATG" and alternative initiation codons
            # 2 = any sense codon
            {"name": "start_pos", "type": "string", "default": '1,2,3'},
            {"name": "genetic_code", "type": "int", "default": 1},
            # Genetic code to use (1-31)
            # see https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi for details
            {"name": "minimal_length", "type": "int", "default": 75},
            # Minimal length of the ORF (nt)
            # Value less than 30 is automatically changed by 30.
            {"name": "strand", "type": "string", "default": "both"},
            # Output ORFs on specified strand only (both|plus|minus)
            {'name': 'update_info', 'type': 'string'},
            {'name': 'main_id', 'type': 'string'}
        ]
        self.add_option(options)
        self.revise_infiles()
        self.tool = self.add_tool("tool_lab.orffinder")
        self.set_options(self._sheet.options())

    def check_options(self):
        if not self.option("input_file").is_set:
            raise OptionError('请输入FASTA文件')
        if self.option("start_codon") not in [0, 1, 2]:
            raise OptionError('开放阅读框参数输入错误')
        if self.option("genetic_code") not in range(1, 32):
            raise OptionError('密码子表参数输入错误')
        return True

    def run(self):
        self.run_tool()
        super(OrffinderWorkflow, self).run()

    def run_tool(self):
        opts = {
            'input_file': self.option('input_file'),
            'start_codon': self.option('start_codon'),
            'start_pos': self.option('start_pos'),
            'genetic_code': self.option('genetic_code'),
            'minimal_length': self.option('minimal_length'),
            'strand': self.option('strand'),
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
        result_dir.add_relpath_rules([
            [".", "", "ORF查找结果文件",0],
            [r'.*.fa', 'xls', 'ORF查找结果文件', 0],
        ])
        super(OrffinderWorkflow, self).end()


class TestFunction(unittest.TestCase):
    '''
    This is test for the workflow. Just run this script to do test.
    '''

    def test(self):
        from mbio.workflows.tool_lab.orffinder import OrffinderWorkflow
        from biocluster.wsheet import Sheet
        import random
        data = {
            "id": "extract_gff_fasta" + str(random.randint(1, 10000)),
            "type": "workflow",
            "name": "tool_lab.extract_gff_fasta",
            "options": dict(
                input_file="/mnt/lustre/users/sanger/sg-users/shicaiping/FASTA_example.fsa",
                start_codon=2,
                start_pos='1',
                genetic_code=1,
                minimal_length=75,
                strand="both",
            )
        }
        wsheet = Sheet(data=data)
        wf =OrffinderWorkflow(wsheet)
        wf.sheet.id = 'orffinder'
        wf.sheet.project_sn = 'orffinder'
        wf.IMPORT_REPORT_DATA = False
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)

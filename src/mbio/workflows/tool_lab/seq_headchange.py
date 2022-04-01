# -*- coding: utf-8 -*-
# __author__ = 'fwy'

import os
from biocluster.workflow import Workflow
import datetime
import unittest
import types
from bson.objectid import ObjectId


class SeqHeadchangeWorkflow(Workflow):
    """
        Used for  change fasta header
        input:
        >seq1
        ATCGTCGATCGA
        >seq1
        TCGAATCGTCATCGTC

        symbol:"seq1"  connect:"_"   position:"back"
        out:
        >S2_seq1
        ATCGTCGATCGA
        >S2_seq1
        TCGAATCGTCATCGTC

        """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(SeqHeadchangeWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "fasta_file", "type": "infile", "format": "ref_rna_v2.common"},  # fasta文件
            # 标识符添加位置  ["back","head"]
            {"name": "position", "type": "string", "default": "before"},
            # 标识符
            {"name": "symbol", "type": "string"},
            # 连接符类型
            {"name": "connect_type", "type": "string"},
            #当标识符选择自定义时客户输入的内容
            {"name": "connect_text", "type": "string"},
            {"name": "main_id", "type": "string",},
            # {"name": "min_len", "type": "int", },
            {'name': 'update_info', 'type': 'string'}
        ]
        self.add_option(options)
        self.revise_infiles()
        self.tool = self.add_tool("tool_lab.seq_headchange")
        self.set_options(self._sheet.options())
        self.connect_dict={"underline":"_","vertical_bar":"|","middle_line":"-"}
        self.connect=self.option("connect_text") if self.option("connect_type")== "custom" else self.connect_dict[self.option("connect_type")]
    def run(self):

        self.fasta_path = self.option("fasta_file").prop["path"]
        self.run_tool()
        super(SeqHeadchangeWorkflow, self).run()

    def run_tool(self):
        opts = {
            'fasta_file': self.fasta_path,
            'position': self.option('position'),
            'symbol': self.option('symbol'),
            'connect': self.connect,
            # 'min_len': self.option('min_len')
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
            [".", "", "fasta序列改ID结果文件",0],
            [r'.*\.cpm2tpm\.xls', 'xls', '定量指标cpm转tpm结果文件', 0],
            [r'.*\.count2tpm\.xls', 'xls', '定量指标count转tpm结果文件', 0],
            [r'.*\.count2cpm\.xls', 'xls', '定量指标count转cpm结果文件', 0],
            [r'.*\.count2fpkm\.xls', 'xls', '定量指标count转fpkm结果文件', 0],
            [r'.*\.fpkm2tpm\.xls', 'xls', '定量指标fpkm转tpm结果文件', 0],
            [r'.*\.cpm2fpkm\.xls', 'xls', '定量指标cpm转fpkm结果文件', 0],
            [r'.*\.fpkm2cpm\.xls', 'xls', '定量指标fpkm转cpm结果文件', 0],
        ])
        super(SeqHeadchangeWorkflow, self).end()


class TestFunction(unittest.TestCase):
    '''
    This is test for the workflow. Just run this script to do test.
    '''

    def test(self):
        from mbio.workflows.tool_lab.seq_headchange import SeqHeadchangeWorkflow
        from biocluster.wsheet import Sheet
        import random
        test_dir = "/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/fungi/Saccharomyces_cerevisiae/Ensemble_release_39/"
        data = {
            "id": "seq_headchange" + str(random.randint(1, 10000)),
            "type": "workflow",
            "name": "tool_lab.seq_headchange",
            "options": dict(
                fasta_file="/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/tool_lab/FAheader/input/ref.fa",
                # table="Standard",
                position="back",
                symbol_type = "underline",
                # symbol_text = "!!!!!",
                # symbol="pppp",
                connect="~"
                # max_len='10000'
            )
        }
        wsheet = Sheet(data=data)
        wf =SeqHeadchangeWorkflow(wsheet)
        wf.sheet.id = 'seq_headchange'
        wf.sheet.project_sn = 'seq_headchange'
        wf.IMPORT_REPORT_DATA = False
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)

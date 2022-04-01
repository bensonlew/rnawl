# -*- coding: utf-8 -*-
# __author__ = 'fwy'

import os
from biocluster.workflow import Workflow
import datetime
import unittest
import types
from bson.objectid import ObjectId


class SeqExtractWorkflow(Workflow):
    """
    Used for seq extract from fasta .

    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(SeqExtractWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "fasta_file", "type": "infile", "format": "ref_rna_v2.common"},  # fasta文件
            {"name": "extract_list", "type": "infile","format":"ref_rna_v2.common" },  #文件格式说明：格式如下，文件中包含五列信息，第一列为染色体名称，第二列为起始位置，第三列为终止位置，第四列为提取序列名称，第五列为链方向性（+正义链/-反义链）。没有表头， 数据与数据之间务必用制表符（tab符）隔开，不能用空格
            {"name": "fasta_seq", "type": "string",},
            {"name": "is_upload", "type": "string","default":"yes" },
            {"name": "main_id", "type": "string"},
            # {"name": "min_len", "type": "int", },
            {'name': 'update_info', 'type': 'string'}
        ]
        self.add_option(options)
        self.revise_infiles()
        self.tool = self.add_tool("tool_lab.seq_extract")
        self.set_options(self._sheet.options())

    def run(self):
        if self.option("is_upload") == "no":
            with open(os.path.join(self.work_dir,"raw_fasta"),"w") as f:
                f.write(self.option("fasta_seq"))
            self.fasta_path =  os.path.join(self.work_dir,"raw_fasta")
        else:
            self.fasta_path = self.option("fasta_file").prop["path"]
        self.run_tool()
        super(SeqExtractWorkflow, self).run()

    def run_tool(self):
        opts = {
            'fasta_file': self.fasta_path,
            'list': self.option('extract_list'),
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
            [".", "", "序列提取结果文件",0],
            # [r'.*\.cpm2tpm\.xls', 'xls', '定量指标cpm转tpm结果文件', 0],
            # [r'.*\.count2tpm\.xls', 'xls', '定量指标count转tpm结果文件', 0],
            # [r'.*\.count2cpm\.xls', 'xls', '定量指标count转cpm结果文件', 0],
            # [r'.*\.count2fpkm\.xls', 'xls', '定量指标count转fpkm结果文件', 0],
            # [r'.*\.fpkm2tpm\.xls', 'xls', '定量指标fpkm转tpm结果文件', 0],
            # [r'.*\.cpm2fpkm\.xls', 'xls', '定量指标cpm转fpkm结果文件', 0],
            # [r'.*\.fpkm2cpm\.xls', 'xls', '定量指标fpkm转cpm结果文件', 0],
        ])
        super(SeqExtractWorkflow, self).end()


class TestFunction(unittest.TestCase):
    '''
    This is test for the workflow. Just run this script to do test.
    '''

    def test(self):
        from mbio.workflows.tool_lab.seq_extract import SeqExtractWorkflow
        from biocluster.wsheet import Sheet
        import random
        test_dir = "/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/fungi/Saccharomyces_cerevisiae/Ensemble_release_39/"
        data = {
            "id": "seq_extract" + str(random.randint(1, 10000)),
            "type": "workflow",
            "name": "tool_lab.seq_extract",
            "options": dict(
                # genome=test_dir + "/" + "dna/Saccharomyces_cerevisiae.dna.toplevel.fa",
                # fasta_seq=">scaffold1\nATTATAAAATGAGGAGATGTAAATTTTAAGGAAAAAAATAAAGCTGTCGAGATTTTCTCGACAGCCTGAAGAGCACCCATCATAAAAGGGTGCTCTTTTTCTTTTCTCTTCCCGATTTGTTCACAGGCTGAATCCTCTCCTCATATGCTGAAAGGGAGATTCAGAAAATATTACAGACTTTTTTAAATAAAGAAGGTGAACGGTCAGTATGTTGAGCGGTTTAACGGTTGCGGTGATCGGGGGAGATGCAAGGCAGCTTGAAATCATTCGCAAGCTGTCACAGCAGCATGCCAAAGTGTTTTTGGTCGGATTTGATCAGCTGGATCATGGGTTTATCGGTGCTGAAAAGCTTAAAATGTCAGAACTTCCATTTGAACAGGTAGACAGTATGATTCTGCCGGTATCAGGTGCAACAGATGAAGGCGTCGTCGCCACAGTTTTCTCAAATGAGCAGGTCGTGCTGGAAGCAGAATATTTAGAAAGAACTCCAGCACATTGTACCTTGTACTCAGGTATTTCTAATACGTACTTAGACAATCTGGCAAAGCAGGTGAACCGGAAGCTTGTGAAGCTGTTTGAGCGCGATGATATTGCCATATATAACTCTATTCCAACAGTTGAAGGGATTATCATGATGGCCATTCAGCAAACGGACTATACGATTCATGGATCACATGTCGCTGTCCTCGGGCTTGGGAGAACAGGGCTCACAATTGCCCGCACAT",
                # is_upload="yes",
                fasta_file="/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/tool_lab/seq_filter/input/example.fasta",
                list="/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/tool_lab/seq_cut/input/seqpos.list"
                # fasta_file="",
                # re_method='c',
                # max_len='10000'
            )
        }
        wsheet = Sheet(data=data)
        wf =SeqExtractWorkflow(wsheet)
        wf.sheet.id = 'seq_extract'
        wf.sheet.project_sn = 'seq_extract'
        wf.IMPORT_REPORT_DATA = False
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)

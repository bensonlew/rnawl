# -*- coding: utf-8 -*-
# __author__ = 'fwy'

import os
from biocluster.workflow import Workflow
import datetime
import unittest
import types
from bson.objectid import ObjectId


class GffposExtractWorkflow(Workflow):
    """
    Used for cds to protein code

    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(GffposExtractWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "gff_file", "type": "infile", "format": "ref_rna_v2.common"},  # 输入文件，可以是gff文件
            {"name": "genelist_file", "type": "infile", "format": "ref_rna_v2.common"},  # 输入文件，gene_id列表
            {"name": "main_id", "type": "string"},
            {'name': 'update_info', 'type': 'string'}
        ]
        self.add_option(options)
        self.revise_infiles()
        self.tool = self.add_tool("tool_lab.gffpos_extract")
        self.set_options(self._sheet.options())

    def run(self):
        self.run_tool()
        super(GffposExtractWorkflow, self).run()


    def run_tool(self):
        opts = {
            'gff_file': self.option('gff_file'),
            'genelist_file': self.option('genelist_file'),
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
             [".", "", "gff位置信息提取结果文件",0],
        #     [r'.*\.cpm2tpm\.xls', 'xls', '定量指标cpm转tpm结果文件', 0],
        #     [r'.*\.count2tpm\.xls', 'xls', '定量指标count转tpm结果文件', 0],
        #     [r'.*\.count2cpm\.xls', 'xls', '定量指标count转cpm结果文件', 0],
        #     [r'.*\.count2fpkm\.xls', 'xls', '定量指标count转fpkm结果文件', 0],
        #     [r'.*\.fpkm2tpm\.xls', 'xls', '定量指标fpkm转tpm结果文件', 0],
        #     [r'.*\.cpm2fpkm\.xls', 'xls', '定量指标cpm转fpkm结果文件', 0],
        #     [r'.*\.fpkm2cpm\.xls', 'xls', '定量指标fpkm转cpm结果文件', 0],
        ])
        super(GffposExtractWorkflow, self).end()


class TestFunction(unittest.TestCase):
    '''
    This is test for the workflow. Just run this script to do test.
    '''

    def test(self):
        from mbio.workflows.tool_lab.gffpos_extract import GffposExtractWorkflow
        from biocluster.wsheet import Sheet
        import random
        test_dir = "/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/fungi/Saccharomyces_cerevisiae/Ensemble_release_39/"
        data = {
            "id": "extract_gff_fasta" + str(random.randint(1, 10000)),
            "type": "workflow",
            "name": "tool_lab.trans_seq",
            "options": dict(
                # genome=test_dir + "/" + "dna/Saccharomyces_cerevisiae.dna.toplevel.fa",
                # fasta_seq=">scaffold1\nATTATAAAATGAGGAGATGTAAATTTTAAGGAAAAAAATAAAGCTGTCGAGATTTTCTCGACAGCCTGAAGAGCACCCATCATAAAAGGGTGCTCTTTTTCTTTTCTCTTCCCGATTTGTTCACAGGCTGAATCCTCTCCTCATATGCTGAAAGGGAGATTCAGAAAATATTACAGACTTTTTTAAATAAAGAAGGTGAACGGTCAGTATGTTGAGCGGTTTAACGGTTGCGGTGATCGGGGGAGATGCAAGGCAGCTTGAAATCATTCGCAAGCTGTCACAGCAGCATGCCAAAGTGTTTTTGGTCGGATTTGATCAGCTGGATCATGGGTTTATCGGTGCTGAAAAGCTTAAAATGTCAGAACTTCCATTTGAACAGGTAGACAGTATGATTCTGCCGGTATCAGGTGCAACAGATGAAGGCGTCGTCGCCACAGTTTTCTCAAATGAGCAGGTCGTGCTGGAAGCAGAATATTTAGAAAGAACTCCAGCACATTGTACCTTGTACTCAGGTATTTCTAATACGTACTTAGACAATCTGGCAAAGCAGGTGAACCGGAAGCTTGTGAAGCTGTTTGAGCGCGATGATATTGCCATATATAACTCTATTCCAACAGTTGAAGGGATTATCATGATGGCCATTCAGCAAACGGACTATACGATTCATGGATCACATGTCGCTGTCCTCGGGCTTGGGAGAACAGGGCTCACAATTGCCCGCACAT",
                # is_upload="yes",
                gff_file="/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/tool_lab_batch2/gtf_bed_gff/data/gff/Homo_sapiens.GRCh38.99.gff3",
                # fasta_file="",
                genelist_file='/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/tool_lab_batch2/gtf_bed_gff/data/gff/testlist',
                # max_len='10000'
            )
        }
        wsheet = Sheet(data=data)
        wf =GffposExtractWorkflow(wsheet)
        wf.sheet.id = 'trans_seq'
        wf.sheet.project_sn = 'trans_seq'
        wf.IMPORT_REPORT_DATA = False
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)

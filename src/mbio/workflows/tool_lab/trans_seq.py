# -*- coding: utf-8 -*-
# __author__ = 'fwy'

import os
from biocluster.workflow import Workflow
import datetime
import unittest
import types
from bson.objectid import ObjectId


class TransSeqWorkflow(Workflow):
    """
    Used for cds to protein code

    detail:
        table:
        OrderedDict([('Standard', '0'),
             ('Standard (with alternative initiation codons)', '1'),
             ('Vertebrate Mitochondrial', '2'),
             ('Yeast Mitochondrial', '3'),
             ('Mold, Protozoan,   Coelenterate Mitochondrial and Mycoplasma/Spiroplasma','4'),
             ('Invertebrate   Mitochondrial', '5'),
             ('Ciliate Macronuclear and  Dasycladacean', '6'),
             ('Echinoderm    Mitochondrial', '9'),
             ('Euplotid Nuclear', '10'),
             ('Bacterial', '11'),
             ('Alternative Yeast Nuclear', '12'),
             ('Ascidian Mitochondrial', '13'),
             ('Flatworm Mitochondrial', '14'),
             ('Blepharisma Macronuclear', '15'),
             ('Chlorophycean Mitochondrial', '16'),
             ('Trematode Mitochondrial', '21'),
             ('Scenedesmus obliquus', '22'),
             ('Thraustochytrium Mitochondrial', '23')])

    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(TransSeqWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "fasta_file", "type": "infile", "format": "ref_rna_v2.common"},  # fasta文件
            {"name": "frame", "type": "int", "default": "1"},
            {"name": "table", "type": "int", "default": "1"},
            #翻译方向 ["5","3","both"]
            {"name": "direction", "type": "string", "default": "5"},
            # 翻译起始位点 ["1","2","3"]
            {"name": "initial_position", "type": "string", "default": "1"},
            {"name": "main_id", "type": "string",}, # 序列转化的方法，默认c，即直接取["c","r","cr"]
            # {"name": "min_len", "type": "int", },
            {'name': 'update_info', 'type': 'string'}
        ]
        self.add_option(options)
        self.revise_infiles()
        self.tool = self.add_tool("tool_lab.trans_seq")
        self.set_options(self._sheet.options())

    def run(self):

        self.fasta_path = self.option("fasta_file").prop["path"]
        self.run_tool()
        super(TransSeqWorkflow, self).run()

    def run_tool(self):
        opts = {
            'fasta_file': self.fasta_path,
            'table': self.option('table'),
            'direction' : self.option("direction"),
            'initial_position':self.option('initial_position')
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
            [".", "", "CDS转蛋白结果文件",0],
            # [r'.*\.cpm2tpm\.xls', 'xls', '定量指标cpm转tpm结果文件', 0],
            # [r'.*\.count2tpm\.xls', 'xls', '定量指标count转tpm结果文件', 0],
            # [r'.*\.count2cpm\.xls', 'xls', '定量指标count转cpm结果文件', 0],
            # [r'.*\.count2fpkm\.xls', 'xls', '定量指标count转fpkm结果文件', 0],
            # [r'.*\.fpkm2tpm\.xls', 'xls', '定量指标fpkm转tpm结果文件', 0],
            # [r'.*\.cpm2fpkm\.xls', 'xls', '定量指标cpm转fpkm结果文件', 0],
            # [r'.*\.fpkm2cpm\.xls', 'xls', '定量指标fpkm转cpm结果文件', 0],
        ])
        super(TransSeqWorkflow, self).end()


class TestFunction(unittest.TestCase):
    '''
    This is test for the workflow. Just run this script to do test.
    '''

    def test(self):
        from mbio.workflows.tool_lab.trans_seq import TransSeqWorkflow
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
                fasta_file="/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/tool_lab/transeq/input/Homo_sapiens.GRCh38.cds.all.fa",
                # fasta_file="",
                table='Standard',
                # max_len='10000'
            )
        }
        wsheet = Sheet(data=data)
        wf =TransSeqWorkflow(wsheet)
        wf.sheet.id = 'trans_seq'
        wf.sheet.project_sn = 'trans_seq'
        wf.IMPORT_REPORT_DATA = False
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)

# -*- coding: utf-8 -*-
# __author__ = 'fwy'

import os
from biocluster.workflow import Workflow
import datetime
import unittest
import types
from bson.objectid import ObjectId


class SeqReverseWorkflow(Workflow):
    """
    Used for fasta file seq reverse .

    detail:
        remethod-r:raw:ATTC new:CTTA
        remethod-c:raw:ATTC new:TAAG
        remethod-cr:raw:ATTC new:GAAT

    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(SeqReverseWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "fasta_file", "type": "infile", "format": "ref_rna_v2.common"},  # fasta文件
            {"name": "fasta_seq", "type": "string",},
            {"name": "is_upload", "type": "string","default":"yes" },
            {"name": "re_method", "type": "string","default":"rc"}, # 序列转化的方法，默认cr，即直接取["c","r","cr"]
            {"name": "main_id", "type": "string"},
            # {"name": "min_len", "type": "int", },
            {'name': 'update_info', 'type': 'string'}
        ]
        self.add_option(options)
        self.revise_infiles()
        self.tool = self.add_tool("tool_lab.seq_reverse")
        self.set_options(self._sheet.options())

    def run(self):
        if self.option("is_upload") == "no":
            with open(os.path.join(self.work_dir,"raw_fasta"),"w") as f:
                f.write(self.option("fasta_seq"))
            self.fasta_path =  os.path.join(self.work_dir,"raw_fasta")
        else:
            self.fasta_path = self.option("fasta_file").prop["path"]
        self.run_tool()
        super(SeqReverseWorkflow, self).run()

    def run_tool(self):
        opts = {
            'fasta_file': self.fasta_path,
            're_method': self.option('re_method'),
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
            [".", "", "序列反转结果文件",0],
        ])
        super(SeqReverseWorkflow, self).end()


class TestFunction(unittest.TestCase):
    '''
    This is test for the workflow. Just run this script to do test.
    '''

    def test(self):
        from mbio.workflows.tool_lab.seq_reverse import SeqReverseWorkflow
        from biocluster.wsheet import Sheet
        import random
        test_dir = "/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/fungi/Saccharomyces_cerevisiae/Ensemble_release_39/"
        data = {
            "id": "extract_gff_fasta" + str(random.randint(1, 10000)),
            "type": "workflow",
            "name": "tool_lab.seq_reverse",
            "options": dict(
                # genome=test_dir + "/" + "dna/Saccharomyces_cerevisiae.dna.toplevel.fa",
                # fasta_seq=">scaffold1\nATTATAAAATGAGGAGATGTAAATTTTAAGGAAAAAAATAAAGCTGTCGAGATTTTCTCGACAGCCTGAAGAGCACCCATCATAAAAGGGTGCTCTTTTTCTTTTCTCTTCCCGATTTGTTCACAGGCTGAATCCTCTCCTCATATGCTGAAAGGGAGATTCAGAAAATATTACAGACTTTTTTAAATAAAGAAGGTGAACGGTCAGTATGTTGAGCGGTTTAACGGTTGCGGTGATCGGGGGAGATGCAAGGCAGCTTGAAATCATTCGCAAGCTGTCACAGCAGCATGCCAAAGTGTTTTTGGTCGGATTTGATCAGCTGGATCATGGGTTTATCGGTGCTGAAAAGCTTAAAATGTCAGAACTTCCATTTGAACAGGTAGACAGTATGATTCTGCCGGTATCAGGTGCAACAGATGAAGGCGTCGTCGCCACAGTTTTCTCAAATGAGCAGGTCGTGCTGGAAGCAGAATATTTAGAAAGAACTCCAGCACATTGTACCTTGTACTCAGGTATTTCTAATACGTACTTAGACAATCTGGCAAAGCAGGTGAACCGGAAGCTTGTGAAGCTGTTTGAGCGCGATGATATTGCCATATATAACTCTATTCCAACAGTTGAAGGGATTATCATGATGGCCATTCAGCAAACGGACTATACGATTCATGGATCACATGTCGCTGTCCTCGGGCTTGGGAGAACAGGGCTCACAATTGCCCGCACAT",
                is_upload="yes",
                fasta_file="/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/tool_lab/seq_filter/input/example.fasta",
                # fasta_file="",
                re_method='c',
                # max_len='10000'
            )
        }
        wsheet = Sheet(data=data)
        wf =SeqReverseWorkflow(wsheet)
        wf.sheet.id = 'extract_gff_fasta'
        wf.sheet.project_sn = 'extract_gff_fasta'
        wf.IMPORT_REPORT_DATA = False
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)

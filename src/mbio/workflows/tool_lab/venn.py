# -*- coding: utf-8 -*-
# __author__ = 'zoujiaxun'

import os
from biocluster.workflow import Workflow
import datetime
import unittest
import types
from bson.objectid import ObjectId


class VennWorkflow(Workflow):
    """
    plot venn

    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(VennWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "venn_file", "type": "infile", "format": "ref_rna_v2.common"},# 用户上传文件
            {'name': 'sep', 'type': 'string'},  # 表格的分隔符
            {'name': 'update_info', 'type': 'string'},
            {'name': "main_id", 'type': 'string'},
        ]
        self.add_option(options)
        self.revise_infiles()
        self.set_options(self._sheet.options())

    def check_options(self):
        # for k, v in self.sheet.options().items():
        #     self.logger.debug('{} = {}'.format(k, v))
        pass

    def run(self):
        self.run_venn()
        super(VennWorkflow, self).run()

    def run_venn(self):
        self.venn = self.add_tool('tool_lab.venn')
        self.venn.set_options({
            'venn_file': self.option('venn_file').prop["path"],
            'sep':self.option('sep')
        })
        self.venn.on('end',self.set_db)
        self.venn.run()


    def end(self):
        super(VennWorkflow, self).end()

    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        venn = self.api.api("tool_lab.venn")
        # add result info
        graph_table = os.path.join(self.venn.output_dir, 'venn_graph.xls')
        venn.add_venn(graph_table, main_id=self.option('main_id'), )
        self.end()

class TestFunction(unittest.TestCase):
    '''
    This is test for the workflow. Just run this script to do test.
    '''

    def test(self):
        from mbio.workflows.tool_lab.venn import VennWorkflow
        from biocluster.wsheet import Sheet
        import random
        # test_dir = "/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/fungi/Saccharomyces_cerevisiae/Ensemble_release_39/"
        data = {
            "id": "venn" + str(random.randint(1, 10000)),
            "type": "workflow",
            "name": "tool_lab.venn",
            "options": dict(
                # genome=test_dir + "/" + "dna/Saccharomyces_cerevisiae.dna.toplevel.fa",
                # fasta_seq=">scaffold1\nATTATAAAATGAGGAGATGTAAATTTTAAGGAAAAAAATAAAGCTGTCGAGATTTTCTCGACAGCCTGAAGAGCACCCATCATAAAAGGGTGCTCTTTTTCTTTTCTCTTCCCGATTTGTTCACAGGCTGAATCCTCTCCTCATATGCTGAAAGGGAGATTCAGAAAATATTACAGACTTTTTTAAATAAAGAAGGTGAACGGTCAGTATGTTGAGCGGTTTAACGGTTGCGGTGATCGGGGGAGATGCAAGGCAGCTTGAAATCATTCGCAAGCTGTCACAGCAGCATGCCAAAGTGTTTTTGGTCGGATTTGATCAGCTGGATCATGGGTTTATCGGTGCTGAAAAGCTTAAAATGTCAGAACTTCCATTTGAACAGGTAGACAGTATGATTCTGCCGGTATCAGGTGCAACAGATGAAGGCGTCGTCGCCACAGTTTTCTCAAATGAGCAGGTCGTGCTGGAAGCAGAATATTTAGAAAGAACTCCAGCACATTGTACCTTGTACTCAGGTATTTCTAATACGTACTTAGACAATCTGGCAAAGCAGGTGAACCGGAAGCTTGTGAAGCTGTTTGAGCGCGATGATATTGCCATATATAACTCTATTCCAACAGTTGAAGGGATTATCATGATGGCCATTCAGCAAACGGACTATACGATTCATGGATCACATGTCGCTGTCCTCGGGCTTGGGAGAACAGGGCTCACAATTGCCCGCACAT",
                venn_file="/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/transcript.count.matrix",
                # fasta_file="/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/tool_lab/seq_filter/input/example.fasta",
                # fasta_file="",
                # target_type="png",
                # max_len='10000'
            )
        }

        wsheet = Sheet(data=data)
        wf =VennWorkflow(wsheet)
        wf.sheet.id = 'venn'
        wf.sheet.project_sn = 'venn'
        wf.IMPORT_REPORT_DATA = False
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)

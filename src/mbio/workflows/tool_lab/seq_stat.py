# -*- coding: utf-8 -*-
# __author__ = 'fwy'

import os
from biocluster.workflow import Workflow
import datetime
import unittest
import types
from bson.objectid import ObjectId


class SeqStatWorkflow(Workflow):
    """
    Used for fasta file seq stat .

    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(SeqStatWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "fasta_file", "type": "infile", "format": "ref_rna_v2.common"},  # fasta文件
            # {"name": "step", "type": "int", "default": 500},  # 统计步长    统计的步长
            {"name": "group_num", "type": "int", "default": 10},  # 统计组数
            # {"name": "min_len", "type": "int", "default": 0},  # 最小统计长度   该长度一下的reads不进行统计
            {"name": "main_id" ,"type": "string"},
            # {"name": "min_len", "type": "int", },
            {'name': 'update_info', 'type': 'string'}
        ]
        self.add_option(options)
        self.revise_infiles()
        self.tool = self.add_tool("tool_lab.seq_stat")
        self.set_options(self._sheet.options())

    def run(self):

        self.run_tool()
        super(SeqStatWorkflow, self).run()

    def run_tool(self):
        opts = {
            'fasta_file': self.option("fasta_file").prop["path"],
            # 'step': self.option('step'),
            'group_num': self.option('group_num')
        }
        self.tool.set_options(opts)
        self.tool.on('end', self.set_db)
        self.tool.run()

    def set_db(self):
        seq_stat = self.api.api("tool_lab.seq_stat")
        seq_stat.add_seq_stat_detail(os.path.join(self.tool.output_dir,"seq_stat.distribution"),os.path.join(self.tool.output_dir,"seq_stat.length"),main_id=self.option("main_id"))
        self.set_output()

    def set_output(self):
        for file in os.listdir(self.tool.output_dir):
            os.link(os.path.join(self.tool.output_dir, file), os.path.join(self.output_dir, file))
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "序列统计结果文件",0],
        ])
        super(SeqStatWorkflow, self).end()


class TestFunction(unittest.TestCase):
    '''
    This is test for the workflow. Just run this script to do test.
    '''

    def test(self):
        from mbio.workflows.tool_lab.seq_stat import SeqStatWorkflow
        from biocluster.wsheet import Sheet
        import random
        test_dir = "/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/fungi/Saccharomyces_cerevisiae/Ensemble_release_39/"
        data = {
            "id": "seq_stat" + str(random.randint(1, 10000)),
            "type": "workflow",
            "name": "tool_lab.seq_stat",
            "options": dict(
                fasta_file="/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/tool_lab/seq_statistic_new/input/example.fasta",
                group_num=10,
                # step=500
            )
        }
        wsheet = Sheet(data=data)
        wf =SeqStatWorkflow(wsheet)
        wf.sheet.id = 'seq_stat'
        wf.sheet.project_sn = 'seq_stat'
        wf.IMPORT_REPORT_DATA = False
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)

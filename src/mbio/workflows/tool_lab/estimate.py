# -*- coding: utf-8 -*-
# __author__ = "chenyanyan, 2016.10.12"
# last_modify by khl 20170504

from biocluster.workflow import Workflow
import os, re, glob
import pandas as pd
import json
from biocluster.core.function import filter_error_info, link, CJsonEncoder
import unittest


class EstimateWorkflow(Workflow):
    """
    差异分析
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(EstimateWorkflow, self).__init__(wsheet_object)
        options = [
            {'name': 'exp_file', 'type': 'infile', 'format': 'ref_rna_v2.common'},
            {'name': 'id', 'type': 'string', 'default': 'GeneSymbol'},
            {'name': 'platform', 'type': 'string', 'default': 'illumina'},
            {'name': 'species', 'type': 'string', 'default': 'Homo_sapiens'},
            {'name': 'sct', 'type': 'string', 'default': 'hierarchy'},
            {'name': 'scm', 'type': 'string', 'default': 'complete'},
            {'name': 'scd', 'type': 'string', 'default': 'correlation'},
            {'name': 'update_info', 'type': 'string'},
            {'name': "main_id", 'type': 'string'},
            {'name': 'sample_str', 'type': 'string', 'default': "All"}
        ]
        self.add_option(options)
        self.revise_infiles()
        self.set_options(self._sheet.options())
        self.estimate = self.add_tool("tool_lab.estimate")


    def run(self):
        self.exp_pre()
        self.run_estimate()
        super(EstimateWorkflow, self).run()

    def exp_pre(self):
        exp = pd.read_table(self.option('exp_file').path, header=0, index_col=0, sep='\t')
        if self.option('sample_str').lower() == 'all' or self.option('sample_str') is None:
            exp_new = exp
        else:
            sample_list = self.option('sample_str').split(';')
            exp_new = exp.loc[:, sample_list]
        self.exp_file_new = os.path.join(self.work_dir, 'count_new.txt')
        exp_new.to_csv(self.exp_file_new, header=True, index=True, sep='\t')
    def run_estimate(self):
        opts = {
            'exp_file': self.exp_file_new,
            'id': self.option('id'),
            'platform': self.option('platform'),
            'species': self.option('species')
        }
        self.estimate.set_options(opts)
        if self.option('sct') == 'no':
            self.estimate.on('end', self.set_db)
        else:
            self.estimate.on('end', self.run_cluster)
        self.estimate.run()

    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        tool_estimate = self.api.api("tool_lab.tool_estimate")
        # add result info
        tool_estimate.add_estimate(self.option('main_id'), self.estimate.option('cluster_matrix').path)
        if self.option('sct') == 'no':
            pass
        else:
            tool_estimate.add_estimate_cluster(self.tool.output_dir, self.option('main_id'))
        self.end()

    def end(self):
        estimate_score_matrix = self.estimate.option('cluster_matrix').path
        os.link(estimate_score_matrix, os.path.join(self.output_dir, 'estimate_score_matrix.txt'))
        result_dir = self.add_upload_dir(self.output_dir)
        super(EstimateWorkflow, self).end()

    def run_cluster(self):
        self.tool = self.add_tool('tool_lab.exp_cluster2gsva')
        options = dict(
            exp=self.estimate.option('cluster_matrix').path,
            sct=self.option('sct'),
            gct='no',
            scm=self.option('scm'),
            gcm='average',
            scd=self.option('scd'),
            gcd='euclidean',
            use_group='no',
        )
        self.tool.set_options(options)
        self.tool.on("end", self.set_db)
        self.tool.run()

class TestFunction(unittest.TestCase):
    '''
    This is test for the workflow. Just run this script to do test.
    '''

    def test(self):
        from mbio.workflows.tool_lab.estimate import EstimateWorkflow
        from biocluster.wsheet import Sheet
        import random
        # test_dir = "/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/fungi/Saccharomyces_cerevisiae/Ensemble_release_39/"
        data = {
            "id": "estimate" + str(random.randint(1, 10000)),
            "type": "workflow",
            "name": "tool_lab.estimate",
            "options": {
                'exp_file': '/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/tool_lab/nomogram/test2medical/count_test.txt',
                'sample_str': 'H1581_1;SNU16_5'
            }
        }

        wsheet = Sheet(data=data)
        wf =EstimateWorkflow(wsheet)
        wf.sheet.id = 'medical_transcriptome'
        wf.sheet.project_sn = 'medical_transcriptome'
        wf.IMPORT_REPORT_DATA = False
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
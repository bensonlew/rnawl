    # -*- coding: utf-8 -*-
# __author__ 'gudeqing,qinjincheng'

from mbio.api.database.ref_rna_v2.api_base import ApiBase
from biocluster.api.database.base import report_check
import re
from collections import defaultdict
import datetime
import unittest
import pandas as pd

class SeqDetail(ApiBase):
    def __init__(self, bind_object):
        super(SeqDetail, self).__init__(bind_object)

    @report_check
    def add_seq_stat(self, gene_stat, trans_stat):
        # task_id = self.bind_object.sheet.id
        # project_sn = self.bind_object.sheet.project_sn
        task_id = self.bind_object.sheet.id
        project_sn = self.bind_object.sheet.project_sn
        time_now = datetime.datetime.now()
        created_ts = time_now.strftime('%Y-%m-%d %H:%M:%S')

        #导入基因序列统计
        name = 'Seq_gene_stat_{}'.format(time_now.strftime('%Y%m%d_%H%M%S'))
        main_info = dict([
            ('task_id', task_id),
            ('project_sn', project_sn),
            ('desc', 'seq_gene_stat'),
            ('created_ts', created_ts),
            ('status', 'start'),
            ('level', 'G'),
            ('name', name),
        ])
        gene_stat_main_id = self.create_db_table('sg_seq_stat', [main_info])
        gene_seq_stat_df = pd.read_table(gene_stat)
        gene_seq_stat_df["seq_stat_id"] = gene_stat_main_id
        gene_seq_stats =gene_seq_stat_df.to_dict("records")
        self.create_db_table('sg_seq_stat_detail', gene_seq_stats)
        self.update_db_record('sg_seq_stat', gene_stat_main_id, status='end', main_id=gene_stat_main_id)

        # 导入转录本序列统计
        name = 'Seq_trans_stat_{}'.format(time_now.strftime('%Y%m%d_%H%M%S'))
        main_info = dict([
            ('task_id', task_id),
            ('project_sn', project_sn),
            ('desc', 'seq_trans_stat'),
            ('created_ts', created_ts),
            ('status', 'start'),
            ('level', 'T'),
            ('name', name),
        ])
        trans_stat_main_id = self.create_db_table('sg_seq_stat', [main_info])
        trans_seq_stat_df = pd.read_table(trans_stat)
        trans_seq_stat_df["seq_stat_id"] = trans_stat_main_id
        trans_seq_stats = trans_seq_stat_df.to_dict("records")
        self.create_db_table('sg_seq_stat_detail', trans_seq_stats)
        self.update_db_record('sg_seq_stat', trans_stat_main_id, status='end', main_id=trans_stat_main_id)

class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """

    def test_gene_detail(self):
        from mbio.workflows.medical_transcriptome.medical_transcriptome_test_api import MedicalTranscriptomeTestApiWorkflow
        from biocluster.wsheet import Sheet
        import random

        data = {
            # "id": "denovo_rna_v2" + str(random.randint(1,10000)),
            "id": "6sieeiseu0ubhjnsm1o6302l5r",
            "project_sn": "45psjlpdvgn6qulm2he6oepva2",
            "type": "workflow",
            "name": "medical_transcriptome.medical_transcriptome_test_api",
            "options": {
            },
        }
        wsheet = Sheet(data=data)
        wf = MedicalTranscriptomeTestApiWorkflow(wsheet)
        wf.IMPORT_REPORT_DATA = True
        wf.IMPORT_REPORT_AFTER_END = False
        wf.test_api = wf.api.api("medical_transcriptome.seq_detail")
        gene_stat="/mnt/ilustre/users/sanger-dev/workspace/20201217/Single_rmats_8392_8879/Detail/output/detail/gene_stat"
        trans_stat="/mnt/ilustre/users/sanger-dev/workspace/20201217/Single_rmats_8392_8879/Detail/output/detail/trans_stat"
        wf.test_api.add_seq_stat(gene_stat, trans_stat)


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTest(TestFunction('test_gene_detail'))
    unittest.TextTestRunner(verbosity=2).run(suite)

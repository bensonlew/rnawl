# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'

import datetime
import os
import pickle as pk
import unittest
import pandas as pd
from api_base import ApiBase


class GeneDetail(ApiBase):
    def __init__(self, bind_object):
        super(GeneDetail, self).__init__(bind_object)

    def add_gene_detail(self, rna_types, result_dir):

        print ", ".join([rna_types, result_dir])
        """
        """
        # create main table

        for rna_type in rna_types.split(","):
            if os.path.exists(os.path.join(result_dir, rna_type)):
                pass
            else:
                raise Exception("{} not found in ".format(os.path.join(result_dir, rna_type)))

        try:
            task_id = self.bind_object.sheet.id
            project_sn = self.bind_object.sheet.project_sn
        except Exception:
            project_sn = task_id = 'do_test'

        create_time = datetime.datetime.now()
        main_info = dict([
            ('task_id', task_id),
            ('project_sn', project_sn),
            ('desc', 'gene_detail_information'),
            ('created_ts', create_time.strftime('%Y-%m-%d %H:%M:%S')),
            ('status', 'start'),
            ('name', 'Gene_detail_' + create_time.strftime("%Y%m%d_%H%M%S")),
            ('refrna_seqdb', result_dir),
            ('version', 'v1.2')
        ])
        main_id = self.create_db_table('genes', [main_info])

        gene_detail_list = list()
        for rna_type in rna_types.split(","):
            if rna_type == "mrna":
                g_detail = os.path.join(result_dir, rna_type, "gene_detail.pk")
                t_detail = os.path.join(result_dir, rna_type, "transcript_detail.pk")
                with open(g_detail) as f:
                    records = pk.load(f)
                    for record in records:
                        record['seq_id'] = record['gene_id']
                        record['level'] = "G"
                        record['category'] = "mRNA"
                        record['genes_id'] = main_id
                        gene_detail_list.append(record)
                with open(t_detail) as f:
                    records = pk.load(f)
                    for record in records:
                        record['seq_id'] = record['transcript_id']
                        record['level'] = "T"
                        record['category'] = "mRNA"
                        record['genes_id'] = main_id
                        gene_detail_list.append(record)
            else:
                t_detail = os.path.join(result_dir, rna_type, "transcript_detail.pk")
                with open(t_detail) as f:
                    records = pk.load(f)
                    for record in records:
                        record['seq_id'] = record['transcript_id']
                        record['level'] = "T"
                        record['category'] = rna_type.replace("rna", "RNA")
                        record['genes_id'] = main_id
                        gene_detail_list.append(record)

        self.create_db_table('genes_detail', gene_detail_list)
        self.update_db_record('genes', main_id, status="end", main_id=main_id)

    def add_seq_stat(self, rna_types, result_dir):

        print ", ".join([rna_types, result_dir])
        """
        """
        # create main table

        for rna_type in rna_types.split(","):
            if os.path.exists(os.path.join(result_dir, rna_type)):
                pass
            else:
                raise Exception("{} not found in ".format(os.path.join(result_dir, rna_type)))

        try:
            task_id = self.bind_object.sheet.id
            project_sn = self.bind_object.sheet.project_sn
        except Exception:
            project_sn = task_id = 'do_test'

        time_now = datetime.datetime.now()
        created_ts = time_now.strftime('%Y-%m-%d %H:%M:%S')
        # 导入基因序列统计
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
        gene_stat_main_id = self.create_db_table('seq_stat', [main_info])
        gene_stat = os.path.join(result_dir, "mrna", "seqdownload_gene_stat")
        gene_seq_stat_df = pd.read_table(gene_stat)
        gene_seq_stat_df["seq_stat_id"] = gene_stat_main_id
        gene_seq_stats = gene_seq_stat_df.to_dict("records")
        self.create_db_table('seq_stat_detail', gene_seq_stats)
        self.update_db_record('seq_stat', gene_stat_main_id, status='end', main_id=gene_stat_main_id)
        self.bind_object.logger.info("导入seq_stat成功!")
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
        trans_stat_main_id = self.create_db_table('seq_stat', [main_info])
        trans_stat = os.path.join(result_dir, "mrna", "seqdownload_tran_stat")
        trans_seq_stat_df = pd.read_table(trans_stat)
        trans_seq_stat_df["seq_stat_id"] = trans_stat_main_id
        trans_seq_stats = trans_seq_stat_df.to_dict("records")
        self.create_db_table('seq_stat_detail', trans_seq_stats)
        self.update_db_record('seq_stat', trans_stat_main_id, status='end', main_id=trans_stat_main_id)
        self.bind_object.logger.info("导入seq_stat成功!")


class TestFunction(unittest.TestCase):
    """
    测试导表函数
    """

    def test_mongo(test):
        from mbio.workflows.whole_transcriptome.whole_transcriptome_test_api import WholeTranscriptomeTestApiWorkflow
        from biocluster.wsheet import Sheet

        data = {
            "id": "whole_transcriptome",
            "project_sn": "whole_transcriptome",
            "type": "workflow",
            "name": "whole_transcriptome.whole_transcriptome_test_api",
            "options": {
            },
        }
        wsheet = Sheet(data=data)
        wf = WholeTranscriptomeTestApiWorkflow(wsheet)

        wf.IMPORT_REPORT_DATA = True
        wf.IMPORT_REPORT_AFTER_END = False
        wf.test_api = wf.api.api("whole_transcriptome.gene_detail")

        test_dir = "/mnt/ilustre/users/sanger-dev/workspace/20191119/Single_whole_5573_9731/GeneDetail/output"
        wf.test_api.add_gene_detail("mrna,mirna,circrna,lncrna", test_dir)


if __name__ == '__main__':
    unittest.main()

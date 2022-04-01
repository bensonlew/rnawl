#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@time    : 2019/2/28 17:52
@file    : lnc_identify.py
"""
# -*- coding: utf-8 -*-
# __author__ = 'wangzhaoyue, shicaiping, qinjincheng'
import csv
import datetime
import json
import os
import unittest

import pandas as pd
from biocluster.api.database.base import report_check
from bson.objectid import ObjectId
from bson.son import SON

from mbio.api.database.lnc_rna.api_base import ApiBase


class LncIdentify(ApiBase):
    def __init__(self, bind_object):
        super(LncIdentify, self).__init__(bind_object)
        self._project_type = 'lnc_rna'

    def create_main_table(self, table_name, content_dict_list):
        main_id = self.create_db_table(table_name, content_dict_list)
        self.update_db_record(table_name, record_id=main_id, main_id=ObjectId(main_id))
        return main_id

    @report_check
    def new_lncrna_predict(self, predictions_detail_path, predictions_stat_path, params=None, tools=None, lnc_annot_path = None):
        """新 lncRNA 预测详情

        :param predictions_dir: 预测结果目录 [lnc_rna.new_lncrna_predict module output目录]
        :param params: 预测参数
        :param tools:
        :return:
        """
        self.bind_object.logger.info('start creating main table in sg_lncrna_new_predict')
        project_sn = self.bind_object.sheet.project_sn
        task_id = self.bind_object.sheet.id
        now_datatime = datetime.datetime.now()
        name = 'LncPredict_' + now_datatime.strftime('%Y%m%d_%H%M%S')
        create_time = now_datatime.strftime('%Y-%m-%d %H:%M:%S')
        insert_data = {
            'project_sn': project_sn,
            'task_id': task_id,
            'name': name,
            # json.dumps(params, sort_keys=True, separators=(',', ':')) separators 为默认值
            'params': params if isinstance(params, str) else json.dumps(params, sort_keys=True),
            'status': 'start',
            'desc': 'new lncrna prediction\'s main table',
            'created_ts': create_time,
        }
        if tools is not None:
            insert_data['tools'] = tools
        main_id = self.create_main_table('sg_new_lncrna_predict', [insert_data])

        self.bind_object.logger.info('succeed in creating main table in sg_new_lncrna_predict')

        annot_ids = set()
        if lnc_annot_path:
            self.bind_object.logger.info('start creating sg_lncrna_new_predict_annot table')
            df = pd.read_table(lnc_annot_path, sep='\t')
            annot_ids = set(df['query_name'])
            data_list = df.to_dict('records')
            self.create_db_table('sg_new_lncrna_predict_annot', data_list, tag_dict={'predict_id': main_id})
            self.bind_object.logger.info('succeed in creating sg_lncrna_new_predict_annot table')


        self.bind_object.logger.info('start creating sg_lncrna_new_predict_detail table')
        df = pd.read_table(predictions_detail_path, sep='\t')
        df.fillna({'gene_description': '', 'gene_name': '', 'gene_id': ''}, inplace=True)
        if lnc_annot_path:
            df['evidence'] = df['transcript_id'].map(lambda x: "yes" if x in annot_ids else "no")
        df = df.round(5)
        data_list = df.to_dict('records')

        self.create_db_table('sg_new_lncrna_predict_detail', data_list, tag_dict={'predict_id': main_id})
        self.bind_object.logger.info('succeed in creating sg_lncrna_new_predict_detail table')
        self.update_db_record('sg_new_lncrna_predict', record_id=main_id, status='end')
        self.bind_object.logger.info('succeed in updating sg_new_lncrna_predict table status')

        self.bind_object.logger.info('start creating main table in sg_new_lncrna_predict_stat')
        insert_data = {
            'project_sn': project_sn,
            'task_id': task_id,
            'name': name,
            'status': 'start',
            'desc': 'new lncrna prediction statistics\'s main table',
            'params': params if isinstance(params, str) else json.dumps(params),
            'created_ts': create_time,
            'tools': tools
        }

        main_id = self.create_main_table('sg_new_lncrna_predict_stat', [insert_data])
        self.bind_object.logger.info('succeed in creating main table in new_lncrna_predict_stat')

        self.bind_object.logger.info('start creating sg_new_lncrna_predict_stat_detail table')

        with open(predictions_stat_path) as in_handler:
            insert_data = [v for v in json.load(in_handler).values()]
        self.create_db_table('sg_new_lncrna_predict_stat_detail', insert_data, tag_dict={'predict_stat_id': main_id})
        self.bind_object.logger.info('succeed in creating sg_new_lncrna_predict_stat_detail table')
        self.update_db_record('sg_new_lncrna_predict_stat', record_id=main_id, status='end')
        self.bind_object.logger.info('succeed in updating sg_new_lncrna_predict_stat table status')

    @report_check
    def known_lncrna_info(self, predictions_detail_path):
        """ 新 lncRNA预测统计 导表

        :param predictions_dir: 预测结果目录 [lnc_rna.new_lncrna_predict module output目录]
        :param params: 预测参数[ 此处忽略 ]
        :return:
        """
        self.bind_object.logger.info('start creating main table in sg_known_lncrna_identify')
        project_sn = self.bind_object.sheet.project_sn
        task_id = self.bind_object.sheet.id
        now_datatime = datetime.datetime.now()
        name = 'KnownLncStat_' + now_datatime.strftime('%Y%m%d_%H%M%S')
        create_time = now_datatime.strftime('%Y-%m-%d %H:%M:%S')
        insert_data = {
            'project_sn': project_sn,
            'task_id': task_id,
            'name': name,
            'status': 'start',
            'desc': 'known lncrna identifying\'s main table',
            'params': 'known lncrna',
            'created_ts': create_time,
        }
        main_id = self.create_main_table('sg_known_lncrna_identify', [insert_data])
        self.bind_object.logger.info('succeed in creating main table in new_lncrna_predict_stat')

        self.bind_object.logger.info('start creating sg_known_lncrna_identify_detail table')

        data_df = pd.read_table(predictions_detail_path, sep='\t', header=0)
        data_df.fillna({'gene_description': '', 'gene_name': '', 'gene_id': ''}, inplace=True)
        detail_data = data_df.to_dict('records')

        self.create_db_table('sg_known_lncrna_identify_detail', detail_data, tag_dict={'identify_id': main_id})
        self.bind_object.logger.info('succeed in creating sg_known_lncrna_identify_detail table')
        self.update_db_record('sg_known_lncrna_identify', record_id=main_id, status='end')
        self.bind_object.logger.info('succeed in updating sg_known_lncrna_identify table status')

    @report_check
    def lncrna_stat(self, lncrna_stat_in_sample, lncrna_stat_in_category):
        """ lncRNA统计 导表

        :param predictions_dir: lncRNA统计结果目录
        :param params: 预测参数[ 此处忽略 ]
        :return:
        """
        self.bind_object.logger.info('start creating main table in sg_lncrna_statistics')
        project_sn = self.bind_object.sheet.project_sn
        task_id = self.bind_object.sheet.id
        now_datatime = datetime.datetime.now()
        name = 'TotalLncStat_' + now_datatime.strftime('%Y%m%d_%H%M%S')
        create_time = now_datatime.strftime('%Y-%m-%d %H:%M:%S')
        insert_data = {
            'project_sn': project_sn,
            'task_id': task_id,
            'name': name,
            'status': 'start',
            'desc': 'known lncrna statistics\'s main table',
            'params': 'known stat',
            'created_ts': create_time,
        }
        main_id = self.create_main_table('sg_lncrna_statistics', [insert_data])
        self.bind_object.logger.info('succeed in creating main table in sg_lncrna_statistics')

        # lncrna 数量在各个样本中统计
        self.bind_object.logger.info('start creating sg_lncrna_statistics_in_samples table')
        data_df = pd.read_table(lncrna_stat_in_sample, sep='\t', header=0)
        detail_data = data_df.to_dict('records')
        self.create_db_table('sg_lncrna_statistics_in_samples', detail_data, tag_dict={'identify_id': main_id})
        self.bind_object.logger.info('succeed in creating sg_lncrna_statistics_in_samples table')

        # lncrna 分类统计
        self.bind_object.logger.info('start creating sg_lncrna_statistics_in_category table')
        data_df = pd.read_table(lncrna_stat_in_category, sep='\t', header=0)
        detail_data = data_df.to_dict('records')
        self.create_db_table('sg_lncrna_statistics_in_category', detail_data, tag_dict={'identify_id': main_id})
        self.bind_object.logger.info('succeed in creating sg_lncrna_statistics_in_category table')

        # 导表状态更新
        self.update_db_record('sg_lncrna_statistics', record_id=main_id, status='end')
        self.bind_object.logger.info('succeed in updating sg_lncrna_statistics table status')


if __name__ == '__main__':
    class TestFunction(unittest.TestCase):
        '''
        This is test for the api. Just run this script to do test.
        '''

        def test(self):
            import random
            from mbio.workflows.lnc_rna.lnc_rna_test_api import LncRnaTestApiWorkflow
            from biocluster.wsheet import Sheet
            data = {
                'id': 'lncrna_predict_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
                'type': 'workflow',
                'name': 'lnc_rna.lnc_rna_test_api',
                'options': {}
            }
            wheet = Sheet(data=data)
            wf = LncRnaTestApiWorkflow(wheet)
            wf.sheet.id = 'lnc_rna'
            wf.sheet.project_sn = 'lnc_rna'
            wf.IMPORT_REPORT_DATA = True
            wf.IMPORT_REPORT_AFTER_DATA = False
            wf.test_api = wf.api.api('lnc_rna.lnc_identify')

            params = dict(
                cpc=True,
                cnci=True,
                cpat=True,
                pfamscan=True,
                # fasta_file=r'/mnt/ilustre/users/sanger-dev/sg-users/zhaozhipeng/lncrna_tools/test_lnc_predict.fa',
                fasta_file=r'/mnt/ilustre/users/isanger/sg-users/qinjincheng/lnc_rna/assemble/'
                           r'method-cufflinks/output/NewTranscripts/new_transcripts.fa',
                gtf_file="/mnt/ilustre/users/isanger/sg-users/qinjincheng/lnc_rna/assemble/"
                         "method-cufflinks/output/NewTranscripts/new_transcripts.gtf",
                identify_num=2,
                transcript_len=200,
                exon_num=2,
                orf_len=300,
                cpc_score=0.5,
                taxonmy='Human',
                cnci_score=0,
                cpat_score=0.5,
                hexamer_dat="/mnt/ilustre/users/isanger/app/bioinfo/lnc_rna/CPAT-1.2.4/dat/Human_Hexamer.tsv",
                logit_model="/mnt/ilustre/users/isanger/app/bioinfo/lnc_rna/CPAT-1.2.4/dat/Human_logitModel.RData",
                new_transcripts_json_file=r"/mnt/ilustre/users/sanger-dev/sg-users/zhaozhipeng/lncrna_tools/new_trans_ids_list.json"
            )

            # wf.test_api.known_lncrna_info(r'/mnt/ilustre/users/sanger-dev/workspace/20190411/Single_new_lncrna_predict1150/KnownLncIdentify/output/known_lncrna_detail.xls')
            # wf.test_api.new_lncrna_predict(
            #     r'/mnt/ilustre/users/sanger-dev/workspace/20190411/Single_MergePredictions_4388/MergePredictions/output/novel_lncrna_predict_detail.xls',
            #     '/mnt/ilustre/users/sanger-dev/workspace/20190411/Single_new_lncrna_predict3512/NewLncrnaPredict/output/novel_lncrna_stat.json',
            #     tools='pfam,cpc,cpat,cnci',
            #     params=params
            # )
            wf.test_api.lncrna_stat('/mnt/ilustre/users/sanger-dev/workspace/20190412/Single_lncrna_stat_5812_7074/LncrnaStat/output/lncrna_stat_in_sample.xls',
                                    '/mnt/ilustre/users/sanger-dev/workspace/20190412/Single_lncrna_stat_5812_7074/LncrnaStat/output/lncrna_stat_in_category.xls')


    unittest.main()

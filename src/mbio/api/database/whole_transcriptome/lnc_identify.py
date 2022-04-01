#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@time    : 2019/2/28 17:52
@file    : lnc_identify.py
"""
# -*- coding: utf-8 -*-
# __author__ = 'wangzhaoyue, shicaiping, qinjincheng'
import datetime
import json
import os
import unittest

import pandas as pd
from biocluster.api.database.base import report_check
from bson.objectid import ObjectId

from mbio.api.database.whole_transcriptome.api_base import ApiBase


class LncIdentify(ApiBase):
    def __init__(self, bind_object):
        super(LncIdentify, self).__init__(bind_object)
        self._project_type = 'whole_transcriptome'

    def create_main_table(self, table_name, content_dict_list):
        main_id = self.create_db_table(table_name, content_dict_list)
        self.update_db_record(table_name, record_id=main_id, main_id=ObjectId(main_id))
        return main_id

    @report_check
    def new_lncrna_predict(self, predictions_detail_path, predictions_stat_path, params=None, tools=None,
                           lnc_annot_path=None):
        """新 lncRNA 预测详情

        :param predictions_dir: 预测结果目录 [whole_transcriptome.new_lncrna_predict module output目录]
        :param params: 预测参数
        :param tools:
        :return:
        """
        self.bind_object.logger.info('start creating main table in lncrna_new_predict')
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
            'version': 'v1'
        }
        if tools is not None:
            insert_data['tools'] = tools
        main_id = self.create_main_table('new_lncrna_predict', [insert_data])

        self.bind_object.logger.info('succeed in creating main table in new_lncrna_predict')

        annot_ids = set()
        if lnc_annot_path:
            self.bind_object.logger.info('start creating lncrna_new_predict_annot table')
            df = pd.read_table(lnc_annot_path, sep='\t')
            annot_ids = set(df['query_name'])
            data_list = df.to_dict('records')
            self.create_db_table('new_lncrna_predict_annot', data_list, tag_dict={'predict_id': main_id})
            self.bind_object.logger.info('succeed in creating lncrna_new_predict_annot table')

        self.bind_object.logger.info('start creating lncrna_new_predict_detail table')
        df = pd.read_table(predictions_detail_path, sep='\t')
        df.fillna({'gene_description': '', 'gene_name': '', 'gene_id': ''}, inplace=True)
        if lnc_annot_path:
            df['evidence'] = df['transcript_id'].map(lambda x: "yes" if x in annot_ids else "no")
        df = df.round(5)
        data_list = df.to_dict('records')

        self.create_db_table('new_lncrna_predict_detail', data_list, tag_dict={'predict_id': main_id})
        self.bind_object.logger.info('succeed in creating lncrna_new_predict_detail table')
        self.update_db_record('new_lncrna_predict', record_id=main_id, status='end')
        self.bind_object.logger.info('succeed in updating new_lncrna_predict table status')

        self.bind_object.logger.info('start creating main table in new_lncrna_predict_stat')
        insert_data = {
            'project_sn': project_sn,
            'task_id': task_id,
            'name': name,
            'status': 'start',
            'desc': 'new lncrna prediction statistics\'s main table',
            'params': params if isinstance(params, str) else json.dumps(params),
            'created_ts': create_time,
            'tools': tools,
            'version': 'v1'
        }

        main_id = self.create_main_table('new_lncrna_predict_stat', [insert_data])
        self.bind_object.logger.info('succeed in creating main table in new_lncrna_predict_stat')

        self.bind_object.logger.info('start creating new_lncrna_predict_stat_detail table')

        with open(predictions_stat_path) as in_handler:
            insert_data = [v for v in json.load(in_handler).values()]
        self.create_db_table('new_lncrna_predict_stat_detail', insert_data, tag_dict={'predict_stat_id': main_id})
        self.bind_object.logger.info('succeed in creating new_lncrna_predict_stat_detail table')
        self.update_db_record('new_lncrna_predict_stat', record_id=main_id, status='end')
        self.bind_object.logger.info('succeed in updating new_lncrna_predict_stat table status')

    @report_check
    def known_lncrna_info(self, predictions_detail_path):
        """ 新 lncRNA预测统计 导表

        :param predictions_dir: 预测结果目录 [whole_transcriptome.new_lncrna_predict module output目录]
        :param params: 预测参数[ 此处忽略 ]
        :return:
        """
        self.bind_object.logger.info('start creating main table in known_lncrna_identify')
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
            'version': 'v1'
        }
        main_id = self.create_main_table('known_lncrna_identify', [insert_data])
        self.bind_object.logger.info('succeed in creating main table in new_lncrna_predict_stat')

        self.bind_object.logger.info('start creating known_lncrna_identify_detail table')

        data_df = pd.read_table(predictions_detail_path, sep='\t', header=0)
        data_df.fillna({'gene_description': '', 'gene_name': '', 'gene_id': ''}, inplace=True)
        detail_data = data_df.to_dict('records')

        self.create_db_table('known_lncrna_identify_detail', detail_data, tag_dict={'identify_id': main_id})
        self.bind_object.logger.info('succeed in creating known_lncrna_identify_detail table')
        self.update_db_record('known_lncrna_identify', record_id=main_id, status='end')
        self.bind_object.logger.info('succeed in updating known_lncrna_identify table status')

    @report_check
    def lncrna_stat(self, lncrna_stat_in_sample, lncrna_stat_in_category, group=None):
        """ lncRNA统计 导表

        :param predictions_dir: lncRNA统计结果目录
        :param params: 预测参数[ 此处忽略 ]
        :return:
        """
        if group:
            sample_list = list()
            with open(group, 'r') as g:
                for line in g.readlines():
                    if line.startswith('#'):
                        continue
                    sample_list.append(line.strip().split('\t')[0])
        self.bind_object.logger.info('start creating main table in lncrna_statistics')
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
            'version': 'v1'
        }
        main_id = self.create_main_table('lncrna_statistics', [insert_data])
        self.bind_object.logger.info('succeed in creating main table in lncrna_statistics')

        # lncrna 数量在各个样本中统计
        self.bind_object.logger.info('start creating lncrna_statistics_in_samples table')
        data_df = pd.read_table(lncrna_stat_in_sample, sep='\t', header=0)
        data_df.sort_values('sample_name', inplace=True)
        detail_data = data_df.to_dict('records')
        if group:
            detail_data.sort(key=lambda x: sample_list.index(x['sample_name']))
        self.create_db_table('lncrna_statistics_in_samples', detail_data, tag_dict={'identify_id': main_id})
        self.bind_object.logger.info('succeed in creating lncrna_statistics_in_samples table')

        # lncrna 分类统计
        self.bind_object.logger.info('start creating lncrna_statistics_in_category table')
        data_df = pd.read_table(lncrna_stat_in_category, sep='\t', header=0)
        detail_data = data_df.to_dict('records')
        self.create_db_table('lncrna_statistics_in_category', detail_data, tag_dict={'identify_id': main_id})
        self.bind_object.logger.info('succeed in creating lncrna_statistics_in_category table')

        # 导表状态更新
        self.update_db_record('lncrna_statistics', record_id=main_id, status='end')
        self.bind_object.logger.info('succeed in updating lncrna_statistics table status')


if __name__ == '__main__':
    class TestFunction(unittest.TestCase):
        '''
        This is test for the api. Just run this script to do test.
        '''

        def test(self):
            import random
            from mbio.workflows.whole_transcriptome.whole_transcriptome_test_api import \
                WholeTranscriptomeTestApiWorkflow
            from biocluster.wsheet import Sheet
            data = {
                'id': 'lncrna_predict_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
                'type': 'workflow',
                'name': 'whole_transcriptome.whole_transcriptome_test_api',
                'options': {}
            }
            wheet = Sheet(data=data)
            wf = WholeTranscriptomeTestApiWorkflow(wheet)
            wf.sheet.id = 'whole_transcriptome_nolncdb'
            wf.sheet.project_sn = 'whole_transcriptome_nolncdb'
            wf.IMPORT_REPORT_DATA = True
            wf.IMPORT_REPORT_AFTER_DATA = False
            wf.test_api = wf.api.api('whole_transcriptome.lnc_identify')

            params = dict(
                cpc=True,
                cnci=True,
                cpat=True,
                pfamscan=True,
                # fasta_file=r'/mnt/ilustre/users/sanger-dev/sg-users/zhaozhipeng/lncrna_tools/test_lnc_predict.fa',
                fasta_file=r'',
                gtf_file="",
                identify_num=2,
                transcript_len=200,
                exon_num=2,
                orf_len=300,
                cpc_score=0.5,
                taxonmy='Human',
                cnci_score=0,
                cpat_score=0.5,
                hexamer_dat="/mnt/ilustre/users/isanger/app/bioinfo/whole_transcriptome/CPAT-1.2.4/dat/Human_Hexamer.tsv",
                logit_model="/mnt/ilustre/users/isanger/app/bioinfo/whole_transcriptome/CPAT-1.2.4/dat/Human_logitModel.RData",
                new_transcripts_json_file=r"/mnt/ilustre/users/sanger-dev/sg-users/zhaozhipeng/lncrna_tools/new_trans_ids_list.json"
            )

            large_gush_dir = "/mnt/ilustre/users/sanger-dev/workspace/20190924/WholeTranscriptome_workflow_1623_1285/LargeGush/output"
            wf.test_api.known_lncrna_info(
                os.path.join(large_gush_dir +
                             "/known_lnc_identify", "known_lncrna_detail.xls"))

            soft_list = [lnc_soft for lnc_soft in ['cpc', 'cnci', 'cpat', 'pfamscan'] if params[lnc_soft]]

            lnc_annot_path = "/mnt/ilustre/users/sanger-dev/workspace/20191008/Single_lnc_annot2144/LncrnaAnnot/output/novel_lncrna_vs_lncrna.xls"
            # lnc_annot_path = None

            wf.test_api.new_lncrna_predict(
                os.path.join(large_gush_dir +
                             "/filter_by_express", "filtered_lncnovel/novel_lncrna_predict_detail.xls"),
                os.path.join(large_gush_dir +
                             "/filter_by_express", "filtered_lncnovel/novel_lncrna_stat.json"),
                tools=",".join(soft_list).replace('pfamscan', 'pfam'),
                params=params,
                lnc_annot_path=lnc_annot_path
            )

            wf.test_api.lncrna_stat(
                os.path.join(large_gush_dir +
                             "/lncrna_stat", 'lncrna_stat_in_sample.xls'),
                os.path.join(large_gush_dir +
                             "/lncrna_stat", 'lncrna_stat_in_category.xls')
            )


    unittest.main()

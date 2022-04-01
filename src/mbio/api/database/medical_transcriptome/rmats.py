# -*- coding: utf-8 -*-
# __author__ = 'shicaiping,qinjincheng'

import datetime
import json
import os
import subprocess
import glob
import pandas as pd
import unittest
from biocluster.api.database.base import report_check
from bson.objectid import ObjectId

from mbio.api.database.medical_transcriptome.api_base import ApiBase



class Rmats(ApiBase):
    def __init__(self, bind_object):
        super(Rmats, self).__init__(bind_object)
        self._RMATS_DETAIL_TABLE_HEAD = \
            ["event_id", "type", "novel", "old", "gene", "gene_symbol", "chr", "strand",
             "exonStart_0base", "exonEnd", "upstreamES", "upstreamEE", "downstreamES",
             "downstreamEE",
             "longExonStart_0base", "longExonEnd", "shortES", "shortEE", "flankingES",
             "flankingEE",
             "riExonStart_0base", "riExonEnd", "1stExonStart_0base", "1stExonEnd",
             "2ndExonStart_0base",
             "2ndExonEnd",
             "diff_JunctionCountOnly", "IJC_SAMPLE_1", "SJC_SAMPLE_1", "IJC_SAMPLE_2",
             "SJC_SAMPLE_2",
             "IncFormLen_JunctionCountOnly", "SkipFormLen_JunctionCountOnly",
             "PValue_JunctionCountOnly",
             "FDR_JunctionCountOnly", "IncLevel1_JunctionCountOnly",
             "IncLevel2_JunctionCountOnly",
             "average_IncLevel1_JunctionCountOnly", "average_IncLevel2_JunctionCountOnly",
             "IncLevelDifference_JunctionCountOnly",
             "increase_inclusion_SAMPLE1_JunctionCountOnly",
             "increase_exclusion_SAMPLE1_JunctionCountOnly",
             "increase_inclusion_SAMPLE2_JunctionCountOnly",
             "increase_exclusion_SAMPLE2_JunctionCountOnly",
             "diff_ReadsOnTargetAndJunctionCounts", "IC_SAMPLE_1", "SC_SAMPLE_1",
             "IC_SAMPLE_2",
             "SC_SAMPLE_2",
             "IncFormLen_ReadsOnTargetAndJunctionCounts",
             "SkipFormLen_ReadsOnTargetAndJunctionCounts",
             "PValue_ReadsOnTargetAndJunctionCounts", "FDR_ReadsOnTargetAndJunctionCounts",
             "IncLevel1_ReadsOnTargetAndJunctionCounts",
             "IncLevel2_ReadsOnTargetAndJunctionCounts",
             "average_IncLevel1_ReadsOnTargetAndJunctionCounts",
             "average_IncLevel2_ReadsOnTargetAndJunctionCounts",
             "IncLevelDifference_ReadsOnTargetAndJunctionCounts",
             "increase_inclusion_SAMPLE1_ReadsOnTargetAndJunctionCounts",
             "increase_exclusion_SAMPLE1_ReadsOnTargetAndJunctionCounts",
             "increase_inclusion_SAMPLE2_ReadsOnTargetAndJunctionCounts",
             "increase_exclusion_SAMPLE2_ReadsOnTargetAndJunctionCounts",
             "diff_JunctionCountOnly_and_ReadsOnTargetAndJunctionCounts",
             "diff_JunctionCountOnly_or_ReadsOnTargetAndJunctionCounts"]

        self._RMATS_DETAIL_MONGO_TABLE_FIELD_DIC = \
            {'event_id': 'event_id', 'type': 'type', 'novel_as': 'novel',
             'gid': 'gene', 'gname': 'gene_symbol', 'chr': 'chr',
             'strand': 'strand', 'es': 'exonStart_0base', 'ee': 'exonEnd',
             'up_es': 'upstreamES', 'up_ee': 'upstreamEE',
             'down_es': 'downstreamES', 'down_ee': 'downstreamEE',
             'les': 'longExonStart_0base', 'lee': 'longExonEnd',
             'ses': 'shortES', 'see': 'shortEE', 'fes': 'flankingES',
             'fee': 'flankingEE', 'ries': 'riExonStart_0base',
             'riee': 'riExonEnd', 'firstes': '1stExonStart_0base',
             'firstee': '1stExonEnd', 'secondes': '2ndExonStart_0base',
             'secondee': '2ndExonEnd', 'diff_jc': 'diff_JunctionCountOnly',
             'ijc_s1': 'IJC_SAMPLE_1', 'sjc_s1': 'SJC_SAMPLE_1',
             'ijc_s2': 'IJC_SAMPLE_2', 'sjc_s2': 'SJC_SAMPLE_2',
             'inclen_jc': 'IncFormLen_JunctionCountOnly',
             'skiplen_jc': 'SkipFormLen_JunctionCountOnly',
             'pvalue_jc': 'PValue_JunctionCountOnly',
             'fdr_jc': 'FDR_JunctionCountOnly',
             'inc1_jc': 'IncLevel1_JunctionCountOnly',
             'inc2_jc': 'IncLevel2_JunctionCountOnly',
             'aver_inc1_jc': 'average_IncLevel1_JunctionCountOnly',
             'aver_inc2_jc': 'average_IncLevel2_JunctionCountOnly',
             'inc_diff_jc': 'IncLevelDifference_JunctionCountOnly',
             'upinc_s1_jc': 'increase_inclusion_SAMPLE1_JunctionCountOnly',
             'upexc_s1_jc': 'increase_exclusion_SAMPLE1_JunctionCountOnly',
             'upinc_s2_jc': 'increase_inclusion_SAMPLE2_JunctionCountOnly',
             'upexc_s2_jc': 'increase_exclusion_SAMPLE2_JunctionCountOnly',
             'diff_all': 'diff_ReadsOnTargetAndJunctionCounts',
             'ic_s1': 'IC_SAMPLE_1', 'sc_s1': 'SC_SAMPLE_1',
             'ic_s2': 'IC_SAMPLE_2', 'sc_s2': 'SC_SAMPLE_2',
             'inclen_all': 'IncFormLen_ReadsOnTargetAndJunctionCounts',
             'skiplen_all': 'SkipFormLen_ReadsOnTargetAndJunctionCounts',
             'pvalue_all': 'PValue_ReadsOnTargetAndJunctionCounts',
             'fdr_all': 'FDR_ReadsOnTargetAndJunctionCounts',
             'inc1_all': 'IncLevel1_ReadsOnTargetAndJunctionCounts',
             'inc2_all': 'IncLevel2_ReadsOnTargetAndJunctionCounts',
             'aver_inc1_all': 'average_IncLevel1_ReadsOnTargetAndJunctionCounts',
             'aver_inc2_all': 'average_IncLevel2_ReadsOnTargetAndJunctionCounts',
             'inc_diff_all': 'IncLevelDifference_ReadsOnTargetAndJunctionCounts',
             'upinc_s1_all': 'increase_inclusion_SAMPLE1_ReadsOnTargetAndJunctionCounts',
             'upexc_s1_all': 'increase_exclusion_SAMPLE1_ReadsOnTargetAndJunctionCounts',
             'upinc_s2_all': 'increase_inclusion_SAMPLE2_ReadsOnTargetAndJunctionCounts',
             'upexc_s2_all': 'increase_exclusion_SAMPLE2_ReadsOnTargetAndJunctionCounts',
             'diff_jc_and_all': 'diff_JunctionCountOnly_and_ReadsOnTargetAndJunctionCounts',
             'diff_jc_or_all': 'diff_JunctionCountOnly_or_ReadsOnTargetAndJunctionCounts'}

        self._RMATS_TYPE_MONGO_TABLE_FIELD_DIC = \
            {'old_a5ss': 'A5SS_old_event_id_set_no',
             'novel_a5ss': 'A5SS_novel_event_id_no', 'total_a5ss': 'A5SS_all_event_id_no',
             'old_se': 'SE_old_event_id_set_no',
             'novel_se': 'SE_novel_event_id_no', 'total_se': 'SE_all_event_id_no',
             'old_mxe': 'MXE_old_event_id_set_no',
             'novel_mxe': 'MXE_novel_event_id_no', 'total_mxe': 'MXE_all_event_id_no',
             'old_ri': 'RI_old_event_id_set_no',
             'novel_ri': 'RI_novel_event_id_no', 'total_ri': 'RI_all_event_id_no',
             'old_a3ss': 'A3SS_old_event_id_set_no',
             'novel_a3ss': 'A3SS_novel_event_id_no', 'total_a3ss': 'A3SS_all_event_id_no',
             'total_total': 'total_as_events_no',
             'novel_total': 'total_as_novel_events_no'}

        self._RMATS_STATS_MONGO_TABLE_FIELD_DIC = \
            {'all_a5ss': 'A5SS_ReadsOnTargetAndJunctionCounts_event_id_set_no',
             'jcandall_a5ss': 'A5SS_JunctionCountOnly_and_ReadsOnTargetAndJunctionCounts_set_no',
             'jc_a5ss': 'A5SS_JunctionCountOnly_event_id_set_no',
             'jcorall_a5ss': 'A5SS_JunctionCountOnly_or_ReadsOnTargetAndJunctionCounts_set_no',
             'all_se': 'SE_ReadsOnTargetAndJunctionCounts_event_id_set_no',
             'jcandall_se': 'SE_JunctionCountOnly_and_ReadsOnTargetAndJunctionCounts_set_no',
             'jc_se': 'SE_JunctionCountOnly_event_id_set_no',
             'jcorall_se': 'SE_JunctionCountOnly_or_ReadsOnTargetAndJunctionCounts_set_no',
             'all_mxe': 'MXE_ReadsOnTargetAndJunctionCounts_event_id_set_no',
             'jcandall_mxe': 'MXE_JunctionCountOnly_and_ReadsOnTargetAndJunctionCounts_set_no',
             'jc_mxe': 'MXE_JunctionCountOnly_event_id_set_no',
             'jcorall_mxe': 'MXE_JunctionCountOnly_or_ReadsOnTargetAndJunctionCounts_set_no',
             'all_ri': 'RI_ReadsOnTargetAndJunctionCounts_event_id_set_no',
             'jcandall_ri': 'RI_JunctionCountOnly_and_ReadsOnTargetAndJunctionCounts_set_no',
             'jc_ri': 'RI_JunctionCountOnly_event_id_set_no',
             'jcorall_ri': 'RI_JunctionCountOnly_or_ReadsOnTargetAndJunctionCounts_set_no',
             'all_a3ss': 'A3SS_ReadsOnTargetAndJunctionCounts_event_id_set_no',
             'jcandall_a3ss': 'A3SS_JunctionCountOnly_and_ReadsOnTargetAndJunctionCounts_set_no',
             'jc_a3ss': 'A3SS_JunctionCountOnly_event_id_set_no',
             'jcorall_a3ss': 'A3SS_JunctionCountOnly_or_ReadsOnTargetAndJunctionCounts_set_no',
             'jc_total': 'total_JunctionCountOnly_event_id_set_no',
             'all_total': 'total_ReadsOnTargetAndJunctionCounts_event_id_set_no'}

        self._RMATS_PSI_MONGO_TABLE_FIELD_DIC = \
            {'a3ss_s1_all_exc': 'A3SS_SAMPLE_1_ReadsOnTargetAndJunctionCounts_exclusion',
             'a3ss_s1_all_inc': 'A3SS_SAMPLE_1_ReadsOnTargetAndJunctionCounts_inclusion',
             'a3ss_s1_all_total': 'A3SS_SAMPLE_1_ReadsOnTargetAndJunctionCounts_total',
             'a3ss_s1_jc_exc': 'A3SS_SAMPLE_1_JunctionCountOnly_exclusion',
             'a3ss_s1_jc_inc': 'A3SS_SAMPLE_1_JunctionCountOnly_inclusion',
             'a3ss_s1_jc_total': 'A3SS_SAMPLE_1_JunctionCountOnly_total',
             'a3ss_s2_all_exc': 'A3SS_SAMPLE_2_ReadsOnTargetAndJunctionCounts_exclusion',
             'a3ss_s2_all_inc': 'A3SS_SAMPLE_2_ReadsOnTargetAndJunctionCounts_inclusion',
             'a3ss_s2_all_total': 'A3SS_SAMPLE_2_ReadsOnTargetAndJunctionCounts_total',
             'a3ss_s2_jc_exc': 'A3SS_SAMPLE_2_JunctionCountOnly_exclusion',
             'a3ss_s2_jc_inc': 'A3SS_SAMPLE_2_JunctionCountOnly_inclusion',
             'a3ss_s2_jc_total': 'A3SS_SAMPLE_2_JunctionCountOnly_total',
             'a5ss_s1_all_exc': 'A5SS_SAMPLE_1_ReadsOnTargetAndJunctionCounts_exclusion',
             'a5ss_s1_all_inc': 'A5SS_SAMPLE_1_ReadsOnTargetAndJunctionCounts_inclusion',
             'a5ss_s1_all_total': 'A5SS_SAMPLE_1_ReadsOnTargetAndJunctionCounts_total',
             'a5ss_s1_jc_exc': 'A5SS_SAMPLE_1_JunctionCountOnly_exclusion',
             'a5ss_s1_jc_inc': 'A5SS_SAMPLE_1_JunctionCountOnly_inclusion',
             'a5ss_s1_jc_total': 'A5SS_SAMPLE_1_JunctionCountOnly_total',
             'a5ss_s2_all_exc': 'A5SS_SAMPLE_2_ReadsOnTargetAndJunctionCounts_exclusion',
             'a5ss_s2_all_inc': 'A5SS_SAMPLE_2_ReadsOnTargetAndJunctionCounts_inclusion',
             'a5ss_s2_all_total': 'A5SS_SAMPLE_2_ReadsOnTargetAndJunctionCounts_total',
             'a5ss_s2_jc_exc': 'A5SS_SAMPLE_2_JunctionCountOnly_exclusion',
             'a5ss_s2_jc_inc': 'A5SS_SAMPLE_2_JunctionCountOnly_inclusion',
             'a5ss_s2_jc_total': 'A5SS_SAMPLE_2_JunctionCountOnly_total',
             'mxe_s1_all_exc': 'MXE_SAMPLE_1_ReadsOnTargetAndJunctionCounts_exclusion',
             'mxe_s1_all_inc': 'MXE_SAMPLE_1_ReadsOnTargetAndJunctionCounts_inclusion',
             'mxe_s1_all_total': 'MXE_SAMPLE_1_ReadsOnTargetAndJunctionCounts_total',
             'mxe_s1_jc_exc': 'MXE_SAMPLE_1_JunctionCountOnly_exclusion',
             'mxe_s1_jc_inc': 'MXE_SAMPLE_1_JunctionCountOnly_inclusion',
             'mxe_s1_jc_total': 'MXE_SAMPLE_1_JunctionCountOnly_total',
             'mxe_s2_all_exc': 'MXE_SAMPLE_2_ReadsOnTargetAndJunctionCounts_exclusion',
             'mxe_s2_all_inc': 'MXE_SAMPLE_2_ReadsOnTargetAndJunctionCounts_inclusion',
             'mxe_s2_all_total': 'MXE_SAMPLE_2_ReadsOnTargetAndJunctionCounts_total',
             'mxe_s2_jc_exc': 'MXE_SAMPLE_2_JunctionCountOnly_exclusion',
             'mxe_s2_jc_inc': 'MXE_SAMPLE_2_JunctionCountOnly_inclusion',
             'mxe_s2_jc_total': 'MXE_SAMPLE_2_JunctionCountOnly_total',
             'ri_s1_all_exc': 'RI_SAMPLE_1_ReadsOnTargetAndJunctionCounts_exclusion',
             'ri_s1_all_inc': 'RI_SAMPLE_1_ReadsOnTargetAndJunctionCounts_inclusion',
             'ri_s1_all_total': 'RI_SAMPLE_1_ReadsOnTargetAndJunctionCounts_total',
             'ri_s1_jc_exc': 'RI_SAMPLE_1_JunctionCountOnly_exclusion',
             'ri_s1_jc_inc': 'RI_SAMPLE_1_JunctionCountOnly_inclusion',
             'ri_s1_jc_total': 'RI_SAMPLE_1_JunctionCountOnly_total',
             'ri_s2_all_exc': 'RI_SAMPLE_2_ReadsOnTargetAndJunctionCounts_exclusion',
             'ri_s2_all_inc': 'RI_SAMPLE_2_ReadsOnTargetAndJunctionCounts_inclusion',
             'ri_s2_all_total': 'RI_SAMPLE_2_ReadsOnTargetAndJunctionCounts_total',
             'ri_s2_jc_exc': 'RI_SAMPLE_2_JunctionCountOnly_exclusion',
             'ri_s2_jc_inc': 'RI_SAMPLE_2_JunctionCountOnly_inclusion',
             'ri_s2_jc_total': 'RI_SAMPLE_2_JunctionCountOnly_total',
             'se_s1_all_exc': 'SE_SAMPLE_1_ReadsOnTargetAndJunctionCounts_exclusion',
             'se_s1_all_inc': 'SE_SAMPLE_1_ReadsOnTargetAndJunctionCounts_inclusion',
             'se_s1_all_total': 'SE_SAMPLE_1_ReadsOnTargetAndJunctionCounts_total',
             'se_s1_jc_exc': 'SE_SAMPLE_1_JunctionCountOnly_exclusion',
             'se_s1_jc_inc': 'SE_SAMPLE_1_JunctionCountOnly_inclusion',
             'se_s1_jc_total': 'SE_SAMPLE_1_JunctionCountOnly_total',
             'se_s2_all_exc': 'SE_SAMPLE_2_ReadsOnTargetAndJunctionCounts_exclusion',
             'se_s2_all_inc': 'SE_SAMPLE_2_ReadsOnTargetAndJunctionCounts_inclusion',
             'se_s2_all_total': 'SE_SAMPLE_2_ReadsOnTargetAndJunctionCounts_total',
             'se_s2_jc_exc': 'SE_SAMPLE_2_JunctionCountOnly_exclusion',
             'se_s2_jc_inc': 'SE_SAMPLE_2_JunctionCountOnly_inclusion',
             'se_s2_jc_total': 'SE_SAMPLE_2_JunctionCountOnly_total',
             'total_s1_all_exc': 'SAMPLE_1_ReadsOnTargetAndJunctionCounts_exclusion_total',
             'total_s1_all_inc': 'SAMPLE_1_ReadsOnTargetAndJunctionCounts_inclusion_total',
             'total_s1_all_total': 'SAMPLE_1_ReadsOnTargetAndJunctionCounts_total',
             'total_s1_jc_exc': 'SAMPLE_1_JunctionCountOnly_exclusion_total',
             'total_s1_jc_inc': 'SAMPLE_1_JunctionCountOnly_inclusion_total',
             'total_s1_jc_total': 'SAMPLE_1_JunctionCountOnly_total',
             'total_s2_all_exc': 'SAMPLE_2_ReadsOnTargetAndJunctionCounts_exclusion_total',
             'total_s2_all_inc': 'SAMPLE_2_ReadsOnTargetAndJunctionCounts_inclusion_total',
             'total_s2_all_total': 'SAMPLE_2_ReadsOnTargetAndJunctionCounts_total',
             'total_s2_jc_exc': 'SAMPLE_2_JunctionCountOnly_exclusion_total',
             'total_s2_jc_inc': 'SAMPLE_2_JunctionCountOnly_inclusion_total',
             'total_s2_jc_total': 'SAMPLE_2_JunctionCountOnly_total'}

    def _parse_dic_file(self, file):
        return dict([(arr[0], arr[1]) for arr in [line.strip().split('\t') for line in open(file).readlines()]])

    @report_check
    def add_sg_splicing_rmats(self, params, outpath):
        task_id = self.bind_object.sheet.id
        project_sn = self.bind_object.sheet.project_sn
        time_now = datetime.datetime.now()
        created_ts = time_now.strftime('%Y-%m-%d %H:%M:%S')
        ctrl, test = os.path.basename(outpath).split('_vs_')
        compare_plan = '{}|{}'.format(ctrl, test)
        name = 'Splicing_{}_vs_{}_{}'.format(ctrl, test, time_now.strftime('%Y%m%d_%H%M%S'))
        result_dir = os.path.join(self.bind_object.sheet.output, '03Gene_structure_analysis','01AS', os.path.basename(outpath)+"_AS")
        rmats_output = os.path.join(self.bind_object.sheet.output.replace('/workflow_results', '/intermediate_results'),
                                   'AS', os.path.basename(outpath))#add by fwy (为配合产品端结果目录修改)
        # result_dir = outpath
        # rmats_output = outpath
        insert_data = {
            'task_id': task_id,
            'project_sn': project_sn,
            'name': name,
            'desc': 'alternative splicing main table',
            'created_ts': created_ts,
            'params': params,
            'status': 'start',
            'result_dir': result_dir,
            'rmats_output': rmats_output,
            'compare_plan': compare_plan,
            "samples": {
                "s1": test,
                "s2": ctrl
            },
            'version': 'v1.1'
        }
        splicing_id = self.create_db_table('sg_splicing_rmats', [insert_data])
        src_detail_file = os.path.join(outpath, 'all_events_detail_big_table.txt')
        type_stats_file = os.path.join(outpath, 'event_type.file.txt')
        chroms = self.add_sg_splicing_rmats_detail(splicing_id=splicing_id, outpath=outpath)[1]
        type_stats_data = self._parse_dic_file(file=type_stats_file)
        self.add_sg_splicing_rmats_type_stats(splicing_id=splicing_id, src_dic=type_stats_data)
        self.add_sg_splicing_rmats_stats(splicing_id=splicing_id, outpath=outpath)
        self.update_db_record('sg_splicing_rmats', splicing_id, status='end',chroms=chroms, main_id=splicing_id)

    @report_check
    def add_sg_splicing_rmats_for_controller(self, splicing_id, outpath, s3_output):
        '''
        {'sg_splicing_rmats': {'sg_splicing_rmats_detail': None,
                               'sg_splicing_rmats_type_stats': None,
                               'sg_splicing_rmats_stats': ('sg_splicing_rmats_diff_stats',
                                                           'sg_splicing_rmats_psi',
                                                           'sg_splicing_rmats_graph')}}
        '''
        splicing_id = ObjectId(splicing_id)
        src_detail_file = os.path.join(outpath, 'all_events_detail_big_table.txt')
        chroms = self.add_sg_splicing_rmats_detail(splicing_id=splicing_id, outpath=outpath)[1]
        type_stats_data = self._parse_dic_file(file=os.path.join(outpath, 'event_type.file.txt'))
        self.add_sg_splicing_rmats_type_stats(splicing_id=splicing_id, src_dic=type_stats_data)
        stat_id = self.add_sg_splicing_rmats_stats(splicing_id=splicing_id, outpath=outpath)
        self.update_db_record('sg_splicing_rmats', splicing_id, status='end', main_id=splicing_id,chroms=chroms,
                              result_dir=os.path.join(s3_output, os.path.basename(outpath)),
                              rmats_output=os.path.join(s3_output, 'rmats_output', os.path.basename(outpath)))
        return stat_id

    def add_sg_splicing_rmats_detail(self, splicing_id, outpath):
        big_table = os.path.join(outpath, 'all_events_detail_big_table.txt')
        event2txpts = dict()
        for event_type in ['SE', 'RI', 'MXE', 'A5SS', 'A3SS']:
            df = pd.read_table(open(os.path.join(outpath, '{}.transcripts.txt'.format(event_type))))
            for i, row in df.iterrows():
                event2txpts.update({row['event_id']: {row.index[1]: row[1], row.index[2]: row[2]}})
        rmats_detail_file_head_index_dic = dict(zip(self._RMATS_DETAIL_TABLE_HEAD, range(64)))
        subprocess.call('sed -i /^[[:space:]]*$/d {}'.format(big_table), shell=True)
        detail_df = pd.read_table(big_table)
        chorms = list(set(detail_df["chr"]))
        file_content = open(big_table)
        file_content.readline()
        documents = list()
        for n, line in enumerate(open(big_table)):
            arr = line.strip().split('\t')
            if n == 0 or len(arr) != 64:
                continue
            data = {'splicing_id': splicing_id}
            for field in self._RMATS_DETAIL_MONGO_TABLE_FIELD_DIC.keys():
                data[field] = arr[rmats_detail_file_head_index_dic[self._RMATS_DETAIL_MONGO_TABLE_FIELD_DIC[field]]]
            try:
                data['inc_diff_jc'] = float(data['inc_diff_jc'])
                data['inc_diff_all'] = float(data['inc_diff_all'])
            except:
                pass
            data['pvalue_jc'] = float(data['pvalue_jc']) if data['pvalue_jc'] != 'null' else data['pvalue_jc']
            data['pvalue_all'] = float(data['pvalue_all']) if data['pvalue_all'] != 'null' else data['pvalue_all']
            data['fdr_jc'] = float(data['fdr_jc']) if data['fdr_jc'] != 'null' else data['fdr_jc']
            data['fdr_all'] = float(data['fdr_all']) if data['fdr_all'] != 'null' else data['fdr_all']
            data['no_diff'] = 'no' if data['diff_jc_or_all'] == 'yes' else 'yes'
            data.update(event2txpts[data['event_id']])
            documents.append(data)
        else:

            collection = self.db['sg_splicing_rmats_detail']
            if documents:
                collection.insert_many(documents)
                return True,chorms

    def add_sg_splicing_rmats_type_stats(self, splicing_id, src_dic):
        dct = dict()
        for field in self._RMATS_TYPE_MONGO_TABLE_FIELD_DIC.keys():
            dct[field] = int(src_dic[self._RMATS_TYPE_MONGO_TABLE_FIELD_DIC[field]])
        dct['old_total'] = int(dct['total_total']) - int(dct['novel_total'])
        as_index = dict(zip(range(3), ['old', 'novel', 'total']))
        document = {'as_stats': {'a3ss': [0, 0, 0], 'a5ss': [0, 0, 0], 'mxe': [0, 0, 0],
                                 'se': [0, 0, 0], 'ri': [0, 0, 0], 'total': [0, 0, 0]},
                    'stats_pie': {'novel': {'A3SS': 0, 'MXE': 0, 'RI': 0, 'SE': 0, 'A5SS': 0},
                                  'old': {'A3SS': 0, 'MXE': 0, 'RI': 0, 'SE': 0, 'A5SS': 0},
                                  'total': {'A3SS': 0, 'MXE': 0, 'RI': 0, 'SE': 0, 'A5SS': 0}},
                    'splicing_id': splicing_id}
        for type in document['as_stats'].keys():
            for i in range(3):
                document['as_stats'][type][i] = int(dct['{}_{}'.format(as_index[i], type)])
        for kind in document['stats_pie'].keys():
            for type in document['stats_pie'][kind].keys():
                document['stats_pie'][kind][type] = dct['{}_{}'.format(kind, type.lower())]
        collection = self.db['sg_splicing_rmats_type_stats']
        collection.insert_one(document)
        return dct

    def add_sg_splicing_rmats_stats(self, splicing_id, outpath):
        document = self.db['sg_splicing_rmats'].find_one({'_id': splicing_id})
        task_id = document['task_id']
        project_sn = document['project_sn']
        compare_plan = document['compare_plan']
        ctrl, test = compare_plan.split('|')
        group = {ctrl: 's1', test: 's2'}
        time_now = datetime.datetime.now()
        name = 'Splicing_stat_{}_vs_{}_{}'.format(ctrl, test, time_now.strftime('%Y%m%d_%H%M%S'))
        created_ts = time_now.strftime('%Y-%m-%d %H:%M:%S')
        params = json.dumps({
            'task_id': task_id,
            'submit_location': 'splicingrmats_stat',
            'task_type': 2,
            'splicing_id': str(splicing_id),
            'psi': '0',
            'pvalue_fdr': 'fdr',
            'fdr': '0.05'
        }, sort_keys=True, separators=(',', ':'))
        insert_data = {
            'task_id': task_id,
            'project_sn': project_sn,
            'name': name,
            'desc': 'alternative splicing stats main table',
            'created_ts': created_ts,
            'params': params,
            'status': 'start',
            'group': group
        }
        collection = self.db['sg_splicing_rmats_stats']
        stat_id = collection.insert_one(insert_data).inserted_id
        diff_stats_data = self._parse_dic_file(file=os.path.join(outpath, 'event_stats.file.txt'))
        diff_stats_dic = self.add_sg_splicing_rmats_diff_stats(stat_id=stat_id, src_dic=diff_stats_data)
        psi_data = self._parse_dic_file(file=os.path.join(outpath, 'psi_stats.file.txt'))
        psi_dic = self.add_sg_splicing_rmats_psi(stat_id=stat_id, src_dic=psi_data)
        self.add_sg_splicing_rmats_graph(stat_id=stat_id, stats=diff_stats_dic, psi=psi_dic)
        self.update_db_record('sg_splicing_rmats_stats', stat_id, status='end', main_id=stat_id)
        return stat_id

    def add_sg_splicing_rmats_stats_for_controller(self, stat_id, outpath):
        stat_id = ObjectId(stat_id)
        diff_stats_data = self._parse_dic_file(file=os.path.join(outpath, 'event_stats.file.txt'))
        diff_stats_dic = self.add_sg_splicing_rmats_diff_stats(stat_id=stat_id, src_dic=diff_stats_data)
        psi_data = self._parse_dic_file(file=os.path.join(outpath, 'psi_stats.file.txt'))
        psi_dic = self.add_sg_splicing_rmats_psi(stat_id=stat_id, src_dic=psi_data)
        self.add_sg_splicing_rmats_graph(stat_id=stat_id, stats=diff_stats_dic, psi=psi_dic)
        self.update_db_record('sg_splicing_rmats_stats', stat_id, status='end', main_id=stat_id)

    def add_sg_splicing_rmats_diff_stats(self, stat_id, src_dic):
        dct = dict()
        for field in self._RMATS_STATS_MONGO_TABLE_FIELD_DIC.keys():
            dct[field] = int(src_dic[self._RMATS_STATS_MONGO_TABLE_FIELD_DIC[field]])
        dct['jcandall_total'] = sum([dct['jcandall_mxe'], dct['jcandall_se'], dct['jcandall_ri'],
                                     dct['jcandall_a3ss'], dct['jcandall_a5ss']])
        dct['jcorall_total'] = sum([dct['jcorall_mxe'], dct['jcorall_se'], dct['jcorall_ri'],
                                    dct['jcorall_a3ss'], dct['jcorall_a5ss']])
        diff_index = dict(zip(range(4), ['jc', 'all', 'jcandall', 'jcorall']))
        document = {'diff_stats': {'a3ss': [0, 0, 0, 0], 'a5ss': [0, 0, 0, 0], 'mxe': [0, 0, 0, 0],
                                   'se': [0, 0, 0, 0], 'ri': [0, 0, 0, 0], 'total': [0, 0, 0, 0]},
                    'stat_id': stat_id}
        for type in document['diff_stats'].keys():
            for i in range(4):
                document['diff_stats'][type][i] = int(dct['{}_{}'.format(diff_index[i], type)])
        collection = self.db['sg_splicing_rmats_diff_stats']
        collection.insert_one(document)
        return dct

    def add_sg_splicing_rmats_psi(self, stat_id, src_dic):
        dct = dict()
        for field in self._RMATS_PSI_MONGO_TABLE_FIELD_DIC.keys():
            dct[field] = int(src_dic[self._RMATS_PSI_MONGO_TABLE_FIELD_DIC[field]])
        document = {'s1_jc': {'a3ss': [0, 0, 0], 'a5ss': [0, 0, 0], 'mxe': [0, 0, 0],
                              'se': [0, 0, 0], 'ri': [0, 0, 0], 'total': [0, 0, 0]},
                    's2_jc': {'a3ss': [0, 0, 0], 'a5ss': [0, 0, 0], 'mxe': [0, 0, 0],
                              'se': [0, 0, 0], 'ri': [0, 0, 0], 'total': [0, 0, 0]},
                    's1_all': {'a3ss': [0, 0, 0], 'a5ss': [0, 0, 0], 'mxe': [0, 0, 0],
                               'se': [0, 0, 0], 'ri': [0, 0, 0], 'total': [0, 0, 0]},
                    's2_all': {'a3ss': [0, 0, 0], 'a5ss': [0, 0, 0], 'mxe': [0, 0, 0],
                               'se': [0, 0, 0], 'ri': [0, 0, 0], 'total': [0, 0, 0]},
                    'stat_id': stat_id}
        psi_index = dict(zip(range(3), ['exc', 'inc', 'total']))
        for choice in ['s1_all', 's1_jc', 's2_all', 's2_jc']:
            for type in document[choice].keys():
                for i in range(3):
                    document[choice][type][i] = dct['{}_{}_{}'.format(type, choice, psi_index[i])]
        collection_obj = self.db['sg_splicing_rmats_psi']
        collection_obj.insert_one(document)
        return dct

    def add_sg_splicing_rmats_graph(self, stat_id, stats, psi):
        document = {'psi_bar': {'s1': {'jc': {'exc': {'A3SS': 0, 'MXE': 0, 'RI': 0, 'SE': 0, 'A5SS': 0},
                                              'inc': {'A3SS': 0, 'MXE': 0, 'RI': 0, 'SE': 0, 'A5SS': 0},
                                              'total': {'A3SS': 0, 'MXE': 0, 'RI': 0, 'SE': 0, 'A5SS': 0}},
                                       'all': {'exc': {'A3SS': 0, 'MXE': 0, 'RI': 0, 'SE': 0, 'A5SS': 0},
                                               'inc': {'A3SS': 0, 'MXE': 0, 'RI': 0, 'SE': 0, 'A5SS': 0},
                                               'total': {'A3SS': 0, 'MXE': 0, 'RI': 0, 'SE': 0, 'A5SS': 0}}},
                                's2': {'jc': {'exc': {'A3SS': 0, 'MXE': 0, 'RI': 0, 'SE': 0, 'A5SS': 0},
                                              'inc': {'A3SS': 0, 'MXE': 0, 'RI': 0, 'SE': 0, 'A5SS': 0},
                                              'total': {'A3SS': 0, 'MXE': 0, 'RI': 0, 'SE': 0, 'A5SS': 0}},
                                       'all': {'exc': {'A3SS': 0, 'MXE': 0, 'RI': 0, 'SE': 0, 'A5SS': 0},
                                               'inc': {'A3SS': 0, 'MXE': 0, 'RI': 0, 'SE': 0, 'A5SS': 0},
                                               'total': {'A3SS': 0, 'MXE': 0, 'RI': 0, 'SE': 0, 'A5SS': 0}}}},
                    'diff_pie': {'jc': {'A3SS': 0, 'MXE': 0, 'RI': 0, 'SE': 0, 'A5SS': 0},
                                 'all': {'A3SS': 0, 'MXE': 0, 'RI': 0, 'SE': 0, 'A5SS': 0},
                                 'jcandall': {'A3SS': 0, 'MXE': 0, 'RI': 0, 'SE': 0, 'A5SS': 0},
                                 'jcorall': {'A3SS': 0, 'MXE': 0, 'RI': 0, 'SE': 0, 'A5SS': 0}},
                    'stat_id': stat_id}
        for src in document['diff_pie'].keys():
            for type in document['diff_pie'][src].keys():
                document['diff_pie'][src][type] = stats['{}_{}'.format(src, type.lower())]
        for sample in document['psi_bar'].keys():
            for src in document['psi_bar'][sample].keys():
                for psi_kind in document['psi_bar'][sample][src].keys():
                    for type in document['psi_bar'][sample][src][psi_kind].keys():
                        document['psi_bar'][sample][src][psi_kind][type] = psi[
                            '_'.join([type.lower(), sample, src, psi_kind])]
        collection = self.db['sg_splicing_rmats_graph']
        collection.insert_one(document)

    def run1(self):
        for p in glob.glob(os.path.join("/mnt/ilustre/users/sanger-dev/workspace/20200812/Single_rmats_6555_4128/Rmats/output", '*')):
            if os.path.isdir(p):
                outpath = p
                ctrl, test = os.path.basename(outpath).split('_vs_')
                # group_dict = {ctrl: self.option('group_table').prop['group_dict'][ctrl],
                #               test: self.option('group_table').prop['group_dict'][test]}
                compare_plan = '{}|{}'.format(ctrl, test)
                params = json.dumps({
                    'task_id': "medical_transcriptome",
                    'submit_location': 'splicingrmats',
                    'task_type': 2,
                    'group_id': "test",
                    'group_dict': "test",
                    'control_id': "test",
                    'compare_plan': compare_plan
                }, sort_keys=True, separators=(',', ':'))
                self.add_sg_splicing_rmats(params, outpath)


class TestFunction(unittest.TestCase):
    '''
    This is test for the api. Just run this script to do test.
    '''

    def test(self):
        import random
        from mbio.workflows.medical_transcriptome.medical_transcriptome_test_api import MedicalTranscriptomeTestApiWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'rmats_count_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'workflow',
            'name': 'medical_transcriptome.medical_transcriptome_test_api',
            'options': {}
        }
        wheet = Sheet(data=data)
        wf = MedicalTranscriptomeTestApiWorkflow(wheet)
        wf.sheet.id = 'medical_transcriptome'
        wf.sheet.project_sn = 'medical_transcriptome'
        wf.IMPORT_REPORT_DATA = True
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.test_api = wf.api.api('medical_transcriptome.rmats')
        params = json.dumps({
            'task_id': "medical_transcriptome",
            'submit_location': 'splicingrmats',
            'task_type': 2,
            'group_id': "5f33947d17b2bf70a4af2eaf",
            'group_dict': {"H1": ["H1581_1", "H1581_2", "H1581_3"], "H2": ["H1581_4", "H1581_5", "H1581_6"]},
            'control_id': "5f33947d17b2bf70a4af2eb0",
            'compare_plan': "H1|H2"
        }, sort_keys=True, separators=(',', ':'))
        wf.test_api.add_sg_splicing_rmats(params,
                                       outpath='/mnt/ilustre/users/sanger-dev/workspace/20200812/Single_rmats_6555_4128/Rmats/output/H1_vs_H2'
                                       )
        params = json.dumps({
            'task_id': "medical_transcriptome",
            'submit_location': 'splicingrmats',
            'task_type': 2,
            'group_id': "5f33947d17b2bf70a4af2eaf",
            'group_dict': {"H1": ["H1581_1", "H1581_2", "H1581_3"], "H3": ["H1581_7", "H1581_8", "H1581_9"]},
            'control_id': "5f33947d17b2bf70a4af2eb0",
            'compare_plan': "H1|H3"
        }, sort_keys=True, separators=(',', ':'))
        wf.test_api.add_sg_splicing_rmats(params,
                                       outpath='/mnt/ilustre/users/sanger-dev/workspace/20200812/Single_rmats_6555_4128/Rmats/output/H1_vs_H3'
                                       )
        params = json.dumps({
            'task_id': "medical_transcriptome",
            'submit_location': 'splicingrmats',
            'task_type': 2,
            'group_id': "5f33947d17b2bf70a4af2eaf",
            'group_dict': {"S1": ["SNU16_1", "SNU16_2", "SNU16_3"], "S2": ["SNU16_4", "SNU16_5", "SNU16_6"]},
            'control_id': "5f33947d17b2bf70a4af2eb0",
            'compare_plan': "S1|S2"
        }, sort_keys=True, separators=(',', ':'))
        wf.test_api.add_sg_splicing_rmats(params,
                                       outpath='/mnt/ilustre/users/sanger-dev/workspace/20200812/Single_rmats_6555_4128/Rmats/output/S1_vs_S2'
                                       )

        params = json.dumps({
            'task_id': "medical_transcriptome",
            'submit_location': 'splicingrmats',
            'task_type': 2,
            'group_id': "5f33947d17b2bf70a4af2eaf",
            'group_dict': {"S1": ["SNU16_1", "SNU16_2", "SNU16_3"], "S3": ["SNU16_7", "SNU16_8", "SNU16_9"]},
            'control_id': "5f33947d17b2bf70a4af2eb0",
            'compare_plan': "S1|S3"
        }, sort_keys=True, separators=(',', ':'))
        wf.test_api.add_sg_splicing_rmats(params,
                                       outpath='/mnt/ilustre/users/sanger-dev/workspace/20200812/Single_rmats_6555_4128/Rmats/output/S1_vs_S3')


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
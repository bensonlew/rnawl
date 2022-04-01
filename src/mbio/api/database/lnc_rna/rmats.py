# -*- coding: utf-8 -*-
# __author__ = 'jinlinfang, shicaiping'

from mbio.api.database.lnc_rna.api_base import ApiBase
from biocluster.api.database.base import Base, report_check
from bson.objectid import ObjectId
import datetime
import json
import os
import re
import subprocess
import unittest

class Rmats(ApiBase):
    '''
    last_modify: 2019.05.10
    '''
    def __init__(self, bind_object):
        super(Rmats, self).__init__(bind_object)
        self._project_type = 'lnc_rna'
        self._RMATS_DETAIL_TABLE_HEAD = [
            "event_id", "type", "novel", "old", "gene", "gene_symbol", "chr", "strand",
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
            "diff_JunctionCountOnly_or_ReadsOnTargetAndJunctionCounts"
        ]
        self._RMATS_DETAIL_MONGO_TABLE_FIELD_DIC = {
            'event_id': 'event_id', 'type': 'type', 'novel_as': 'novel',
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
            'diff_jc_or_all': 'diff_JunctionCountOnly_or_ReadsOnTargetAndJunctionCounts'
        }
        self._RMATS_TYPE_MONGO_TABLE_FIELD_DIC = {
            'old_a5ss': 'A5SS_old_event_id_set_no',
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
            'novel_total': 'total_as_novel_events_no'
        }
        self._RMATS_STATS_MONGO_TABLE_FIELD_DIC = {
            'all_a5ss': 'A5SS_ReadsOnTargetAndJunctionCounts_event_id_set_no',
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
            'all_total': 'total_ReadsOnTargetAndJunctionCounts_event_id_set_no',
        }
        self._RMATS_PSI_MONGO_TABLE_FIELD_DIC = {
            'a3ss_s1_all_exc': 'A3SS_SAMPLE_1_ReadsOnTargetAndJunctionCounts_exclusion',
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
            'total_s2_jc_total': 'SAMPLE_2_JunctionCountOnly_total'
        }

    def parse_dic_file(self, dic_file):
        return dict([(arr[0], arr[1]) for arr in [line.strip().split('\t') for line in open(dic_file).readlines()]])

    @report_check
    def add_sg_splicing_rmats(self, splicing_id=None, outpath=None, compare_plan=None, params=None, s3_output=None, gene_type_tsv=None):
        if splicing_id:
            splicing_id = ObjectId(splicing_id)
        else:
            project_sn = self.bind_object.sheet.project_sn
            task_id = self.bind_object.sheet.id
            test_group = compare_plan.split('|')[0]
            ctrl_group = compare_plan.split('|')[1]
            time_now = datetime.datetime.now()
            name = 'Splicing_{}_vs_{}_{}'.format(test_group, ctrl_group, time_now.strftime('%Y%m%d_%H%M%S'))
            insert_data = {
                'project_sn': project_sn,
                'task_id': task_id,
                'name': name,
                'desc': 'alternative splicing main table',
                'created_ts': time_now.strftime('%Y-%m-%d %H:%M:%S'),
                'params': json.dumps(params, sort_keys=True, separators=(',', ':')),
                'status': 'start',
                'result_dir': outpath,
                'compare_plan': compare_plan
            }
            splicing_id = self.create_db_table('sg_splicing_rmats', [insert_data])
            self.bind_object.logger.info('succeed in creating table in sg_splicing_rmats')
        src_detail_file = os.path.join(outpath, 'all_events_detail_big_table.txt')
        self.add_sg_splicing_rmats_detail(splicing_id=splicing_id, src_file=src_detail_file, type_file=gene_type_tsv)
        type_stats_data = self.parse_dic_file(dic_file=os.path.join(outpath, 'event_type.file.txt'))
        self.add_sg_splicing_rmats_type_stats(splicing_id=splicing_id, src_dic=type_stats_data)
        self.add_sg_splicing_rmats_stats(splicing_id=splicing_id, outpath=outpath)
        insert_dict = {'status': 'end', 'main_id': splicing_id}
        if s3_output:
            insert_dict.update({'result_dir': os.path.join(s3_output, os.path.basename(outpath))})
        self.update_db_record('sg_splicing_rmats', splicing_id, insert_dict=insert_dict)

    def add_sg_splicing_rmats_detail(self, splicing_id, src_file, type_file):
        gene_type_dict = dict()
        for line in open(type_file):
            items = line.strip().split('\t')
            if len(items) > 3:
                gene_type_dict[items[0]] = items[2]
        rmats_detail_file_head_index_dic = dict(zip(self._RMATS_DETAIL_TABLE_HEAD, range(64)))
        subprocess.call('sed -i /^[[:space:]]*$/d %s' % (src_file), shell=True)
        file_content = open(src_file)
        file_content.readline()
        data_lst = []
        while 1:
            data = {}
            line = file_content.readline()
            if not line.strip():
                break
            arr = re.split('\t+', line.strip())
            if len(arr) != 64:
                raise Exception('line in rmats big detail file %s is not legal: not 64 columns: %s' % (src_file, line))
            data['splicing_id'] = splicing_id
            for field in self._RMATS_DETAIL_MONGO_TABLE_FIELD_DIC.keys():
                data[field] = arr[rmats_detail_file_head_index_dic[self._RMATS_DETAIL_MONGO_TABLE_FIELD_DIC[field]]]

            try:
                data['inc_diff_jc'] = float(data['inc_diff_jc'])
                data['inc_diff_all'] = float(data['inc_diff_all'])
            except:
                pass
            data["pvalue_jc"] = float(data["pvalue_jc"]) if data["pvalue_jc"] != "null" else data["pvalue_jc"]
            data["pvalue_all"] = float(data["pvalue_all"]) if data["pvalue_all"] != "null" else data["pvalue_all"]
            data["fdr_jc"] = float(data["fdr_jc"]) if data["fdr_jc"] != "null" else data["fdr_jc"]
            data["fdr_all"] = float(data["fdr_all"]) if data["fdr_all"] != "null" else data["fdr_all"]
            data['splicing_id'] = splicing_id
            data['no_diff'] = 'no' if data['diff_jc_or_all'] == 'yes' else 'yes'
            if data['gid'] in gene_type_dict:
                data.update({'rna_type': gene_type_dict[data['gid']]})
            else:
                data.update({'rna_type': 'other'})
            data_lst.append(data)
        collection_obj = self.db['sg_splicing_rmats_detail']
        file_content = open(src_file)
        if len(file_content.readlines()) >= 2:
            collection_obj.insert_many(data_lst)
            self.bind_object.logger.info('succeed in creating tables in sg_splicing_rmats_detail')
        return True

    def add_sg_splicing_rmats_stats(self, splicing_id=None, outpath=None):
        task_id = self.db['sg_splicing_rmats'].find_one({"_id": splicing_id})["task_id"]
        compare_plan = self.db['sg_splicing_rmats'].find_one({"_id": splicing_id})['compare_plan']
        test_group = compare_plan.split('|')[0]
        ctrl_group = compare_plan.split('|')[1]
        group = {test_group: 's1', ctrl_group: 's2'}
        time_now = datetime.datetime.now()
        name = 'Splicing_stat_{}_vs_{}_{}'.format(test_group, ctrl_group, time_now.strftime('%Y%m%d_%H%M%S'))
        params = {
            'task_id': task_id,
            'submit_location': 'splicingrmats_stat',
            'task_type': 2,
            'splicing_id': str(splicing_id),
            'delta_psi': '0',
            'significant_diff': 'fdr',
            'significant_value': '0.05',
        }
        insert_data = {
            'task_id': task_id,
            'name': name,
            'desc': 'alternative splicing statistics main table',
            'created_ts': time_now.strftime('%Y-%m-%d %H:%M:%S'),
            'params': json.dumps(params, sort_keys=True, separators=(',', ':')),
            'status': 'start',
            'group': group
        }
        collection_obj = self.db['sg_splicing_rmats_stats']
        stat_id = collection_obj.insert_one(insert_data).inserted_id
        if stat_id:
            self.bind_object.logger.info('succeed in creating table in sg_splicing_rmats_stats')
        self.update_db_record('sg_splicing_rmats_stats', stat_id, status='end', main_id=stat_id)
        diff_stats_file = os.path.join(outpath, 'event_stats.file.txt')
        psi_stat_file = os.path.join(outpath, 'psi_stats.file.txt')
        psi_data = self.parse_dic_file(dic_file=psi_stat_file)
        diff_stats_data = self.parse_dic_file(dic_file=diff_stats_file)
        diff_stats_dic = self.add_sg_splicing_rmats_diff_stats(stat_id=stat_id, src_dic=diff_stats_data, outpath=outpath)
        psi_dic = self.add_sg_splicing_rmats_psi(stat_id=stat_id, src_dic=psi_data)
        rmats_graph_id = self.add_sg_splicing_rmats_graph(stat_id=stat_id, stats=diff_stats_dic, psi=psi_dic, outpath=outpath)
        if rmats_graph_id:
            self.update_db_record('sg_splicing_rmats_stats', stat_id, status='end')

    def add_sg_splicing_rmats_type_stats(self, splicing_id, src_dic):
        tmp_dic = {}
        for field in self._RMATS_TYPE_MONGO_TABLE_FIELD_DIC.keys():
            tmp_dic[field] = int(src_dic[self._RMATS_TYPE_MONGO_TABLE_FIELD_DIC[field]])
        tmp_dic['old_total'] = int(tmp_dic['total_total']) - int(tmp_dic['novel_total'])
        as_index = dict(zip(range(3), ['old', 'novel', 'total']))
        stats_dic = {
            'as_stats': {'a3ss': [0, 0, 0], 'a5ss': [0, 0, 0], 'mxe': [0, 0, 0], 'se': [0, 0, 0], 'ri': [0, 0, 0],
                         'total': [0, 0, 0]},
            'stats_pie': {'novel': {'A3SS': 0, 'MXE': 0, 'RI': 0, 'SE': 0, 'A5SS': 0, },
                          'old': {'A3SS': 0, 'MXE': 0, 'RI': 0, 'SE': 0, 'A5SS': 0, },
                          'total': {'A3SS': 0, 'MXE': 0, 'RI': 0, 'SE': 0, 'A5SS': 0, }},
            'splicing_id': splicing_id
        }
        for type in stats_dic['as_stats'].keys():
            for i in range(3):
                stats_dic['as_stats'][type][i] = int(tmp_dic[as_index[i] + '_' + type])
        for kind in stats_dic['stats_pie'].keys():
            for type in stats_dic['stats_pie'][kind].keys():
                stats_dic['stats_pie'][kind][type] = tmp_dic[kind + '_' + type.lower()]
        collection_obj = self.db['sg_splicing_rmats_type_stats']
        _id = collection_obj.insert_one(stats_dic).inserted_id
        if _id:
            self.bind_object.logger.info('succeed in creating table in sg_splicing_rmats_type_stats')
        return tmp_dic

    def add_sg_splicing_rmats_diff_stats(self, stat_id, src_dic, outpath=None):
        tmp_dic = {}
        for field in self._RMATS_STATS_MONGO_TABLE_FIELD_DIC.keys():
            tmp_dic[field] = int(src_dic[self._RMATS_STATS_MONGO_TABLE_FIELD_DIC[field]])
        tmp_dic['jcandall_total'] = sum(
            [tmp_dic['jcandall_mxe'], tmp_dic['jcandall_se'], tmp_dic['jcandall_ri'],
             tmp_dic['jcandall_a3ss'], tmp_dic['jcandall_a5ss']])
        tmp_dic['jcorall_total'] = sum(
            [tmp_dic['jcorall_mxe'], tmp_dic['jcorall_se'], tmp_dic['jcorall_ri'],
             tmp_dic['jcorall_a3ss'], tmp_dic['jcorall_a5ss']])
        diff_index = dict(zip(range(4), ['jc', 'all', 'jcandall', 'jcorall']))
        stats_dic = {
            'diff_stats':{'a3ss': [0, 0, 0, 0], 'a5ss': [0, 0, 0, 0], 'mxe': [0, 0, 0, 0], 'se': [0, 0, 0, 0],
                          'ri': [0, 0, 0, 0], 'total': [0, 0, 0, 0]},
            'stat_id': stat_id,
        }
        for type in stats_dic['diff_stats'].keys():
            for i in range(4):
                stats_dic['diff_stats'][type][i] = int(tmp_dic[diff_index[i] + '_' + type])

        # export to excel on 20201221
        splicing_stats = os.path.join(outpath, 'splicing_stats.xls')
        with open(splicing_stats, 'w') as ss:
            ss.write('AS type\tJunctionCountOnly(JC)\tReadsOnTargetAndJunctionCounts(JCEC)\tJC&JCEC\tJC|JCEC\n')
            for key, value in stats_dic['diff_stats'].items():
                ss.write('{}\t{}\n'.format(key, '\t'.join([str(i) for i in value])))

        collection_obj = self.db['sg_splicing_rmats_diff_stats']
        _id = collection_obj.insert_one(stats_dic).inserted_id
        if _id:
            self.bind_object.logger.info('succeed in creating table in sg_splicing_rmats_diff_stats')
        return tmp_dic

    def add_sg_splicing_rmats_psi(self, stat_id, src_dic):
        tmp_dic = {}
        for field in self._RMATS_PSI_MONGO_TABLE_FIELD_DIC.keys():
            tmp_dic[field] = int(src_dic[self._RMATS_PSI_MONGO_TABLE_FIELD_DIC[field]])
        psi_dic = {
            's1_jc':
                {'a3ss': [0, 0, 0], 'a5ss': [0, 0, 0], 'mxe': [0, 0, 0], 'se': [0, 0, 0], 'ri': [0, 0, 0],
                 'total': [0, 0, 0]
                 },
            's2_jc':
                {'a3ss': [0, 0, 0], 'a5ss': [0, 0, 0], 'mxe': [0, 0, 0], 'se': [0, 0, 0], 'ri': [0, 0, 0],
                 'total': [0, 0, 0]
                 },
            's1_all':
                {'a3ss': [0, 0, 0], 'a5ss': [0, 0, 0], 'mxe': [0, 0, 0], 'se': [0, 0, 0], 'ri': [0, 0, 0],
                 'total': [0, 0, 0]
                 },
            's2_all':
                {'a3ss': [0, 0, 0], 'a5ss': [0, 0, 0], 'mxe': [0, 0, 0], 'se': [0, 0, 0], 'ri': [0, 0, 0],
                 'total': [0, 0, 0]
                 },
            'stat_id': stat_id
        }
        psi_index = dict(zip(range(3), ['exc', 'inc', 'total']))
        for choice in ['s1_all', 's1_jc', 's2_all', 's2_jc']:
            for type in psi_dic[choice].keys():
                for i in range(3):
                    psi_dic[choice][type][i] = tmp_dic[type + '_' + choice + '_' + psi_index[i]]

        collection_obj = self.db['sg_splicing_rmats_psi']
        _id = collection_obj.insert_one(psi_dic).inserted_id
        if _id:
            self.bind_object.logger.info('succeed in creating table in sg_splicing_rmats_psi')

        return tmp_dic

    def add_sg_splicing_rmats_graph(self, stat_id, stats, psi, outpath=None):
        graph_dic = {'psi_bar': {'s1': {'jc': {'exc': {'A3SS': 0, 'MXE': 0, 'RI': 0, 'SE': 0, 'A5SS': 0, },
                                               'inc': {'A3SS': 0, 'MXE': 0, 'RI': 0, 'SE': 0, 'A5SS': 0, },
                                               'total': {'A3SS': 0, 'MXE': 0, 'RI': 0, 'SE': 0, 'A5SS': 0, }
                                               }, 'all': {'exc': {'A3SS': 0, 'MXE': 0, 'RI': 0, 'SE': 0, 'A5SS': 0, },
                                                          'inc': {'A3SS': 0, 'MXE': 0, 'RI': 0, 'SE': 0, 'A5SS': 0, },
                                                          'total': {'A3SS': 0, 'MXE': 0, 'RI': 0, 'SE': 0,
                                                                    'A5SS': 0, }}},
                                 's2': {'jc': {'exc': {'A3SS': 0, 'MXE': 0, 'RI': 0, 'SE': 0, 'A5SS': 0, },
                                               'inc': {'A3SS': 0, 'MXE': 0, 'RI': 0, 'SE': 0, 'A5SS': 0, },
                                               'total': {'A3SS': 0, 'MXE': 0, 'RI': 0, 'SE': 0, 'A5SS': 0, }},
                                        'all': {'exc': {'A3SS': 0, 'MXE': 0, 'RI': 0, 'SE': 0, 'A5SS': 0, },
                                                'inc': {'A3SS': 0, 'MXE': 0, 'RI': 0, 'SE': 0, 'A5SS': 0, },
                                                'total': {'A3SS': 0, 'MXE': 0, 'RI': 0, 'SE': 0, 'A5SS': 0, }}}},
                     'diff_pie': {'jc': {'A3SS': 0, 'MXE': 0, 'RI': 0, 'SE': 0, 'A5SS': 0, },
                                  'all': {'A3SS': 0, 'MXE': 0, 'RI': 0, 'SE': 0, 'A5SS': 0, },
                                  'jcandall': {'A3SS': 0, 'MXE': 0, 'RI': 0, 'SE': 0, 'A5SS': 0, },
                                  'jcorall': {'A3SS': 0, 'MXE': 0, 'RI': 0, 'SE': 0, 'A5SS': 0, }},
                     'stat_id': stat_id}

        for src in graph_dic['diff_pie'].keys():
            for type in graph_dic['diff_pie'][src].keys():
                graph_dic['diff_pie'][src][type] = stats[src + '_' + type.lower()]
        for sample in graph_dic['psi_bar'].keys():
            for src in graph_dic['psi_bar'][sample].keys():
                for psi_kind in graph_dic['psi_bar'][sample][src].keys():
                    for type in graph_dic['psi_bar'][sample][src][psi_kind].keys():
                        graph_dic['psi_bar'][sample][src][psi_kind][type] = psi[
                            '_'.join([type.lower(), sample, src, psi_kind])]

        # export to excel on 20201215
        as_type = ['A3SS', 'MXE', 'RI', 'SE', 'A5SS']
        main_collection = self.db['sg_splicing_rmats_stats']
        sample_dict = main_collection.find_one({'_id': stat_id})['group']
        for sample in graph_dic['psi_bar'].keys():
            for k, v in sample_dict.items():    # get sample name
                if sample == v:
                    sample_name = k
            for src in graph_dic['psi_bar'][sample].keys():
                if src == 'all':                  # rename diff method
                    method = 'JCEC'
                else:
                    method = src.upper()
                stats_file = os.path.join(outpath, '{}_splicing_{}_stats.xls'.format(sample_name, method))
                with open(stats_file, 'w') as stats:
                    stats.write('AS type\tExclusion\tInclusion\tTotal events\n')
                    for i in as_type:
                        bar_list = list()
                        for psi_kind in graph_dic['psi_bar'][sample][src].keys():
                            bar_list.append(str(graph_dic['psi_bar'][sample][src][psi_kind].get(i)))
                        stats.write('{}\t{}\n'.format(i, '\t'.join(bar_list)))

        collection_obj = self.db['sg_splicing_rmats_graph']
        _id = collection_obj.insert_one(graph_dic).inserted_id
        if _id:
            self.bind_object.logger.info('succeed in creating table in sg_splicing_rmats_graph')
            return _id

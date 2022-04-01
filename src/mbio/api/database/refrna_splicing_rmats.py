# -*- coding: utf-8 -*-
# __author__ = fiona
# time: 2017/5/9 14:14

import re, os, Bio, argparse, sys, fileinput, urllib2

from pymongo import MongoClient
from bson.objectid import ObjectId
import types
from types import StringTypes
import re, subprocess
import json, time
import pandas as pd
import numpy as np
import datetime, os
from bson.son import SON
from collections import Counter
import glob
from biocluster.api.database.base import Base, report_check
from biocluster.config import Config


class RefrnaSplicingRmats(Base):
    def __init__(self, bind_object):
        super(RefrnaSplicingRmats, self).__init__(bind_object)
        self._project_type = 'ref_rna'
        #self._db_name = Config().MONGODB + '_ref_rna'
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
        
        self._RMATS_STATS_MONGO_TABLE_FIELD_DIC = \
            {
                'all_a5ss': 'A5SS_ReadsOnTargetAndJunctionCounts_event_id_set_no',
                'jcandall_a5ss': 'A5SS_JunctionCountOnly_and_ReadsOnTargetAndJunctionCounts_set_no',
                'old_a5ss': 'A5SS_old_event_id_set_no', 'jc_a5ss': 'A5SS_JunctionCountOnly_event_id_set_no',
                'novel_a5ss': 'A5SS_novel_event_id_no', 'total_a5ss': 'A5SS_all_event_id_no',
                'jcorall_a5ss': 'A5SS_JunctionCountOnly_or_ReadsOnTargetAndJunctionCounts_set_no',
                'all_se': 'SE_ReadsOnTargetAndJunctionCounts_event_id_set_no',
                'jcandall_se': 'SE_JunctionCountOnly_and_ReadsOnTargetAndJunctionCounts_set_no',
                'old_se': 'SE_old_event_id_set_no', 'jc_se': 'SE_JunctionCountOnly_event_id_set_no',
                'novel_se': 'SE_novel_event_id_no', 'total_se': 'SE_all_event_id_no',
                'jcorall_se': 'SE_JunctionCountOnly_or_ReadsOnTargetAndJunctionCounts_set_no',
                'all_mxe': 'MXE_ReadsOnTargetAndJunctionCounts_event_id_set_no',
                'jcandall_mxe': 'MXE_JunctionCountOnly_and_ReadsOnTargetAndJunctionCounts_set_no',
                'old_mxe': 'MXE_old_event_id_set_no', 'jc_mxe': 'MXE_JunctionCountOnly_event_id_set_no',
                'novel_mxe': 'MXE_novel_event_id_no', 'total_mxe': 'MXE_all_event_id_no',
                'jcorall_mxe': 'MXE_JunctionCountOnly_or_ReadsOnTargetAndJunctionCounts_set_no',
                'all_ri': 'RI_ReadsOnTargetAndJunctionCounts_event_id_set_no',
                'jcandall_ri': 'RI_JunctionCountOnly_and_ReadsOnTargetAndJunctionCounts_set_no',
                'old_ri': 'RI_old_event_id_set_no', 'jc_ri': 'RI_JunctionCountOnly_event_id_set_no',
                'novel_ri': 'RI_novel_event_id_no', 'total_ri': 'RI_all_event_id_no',
                'jcorall_ri': 'RI_JunctionCountOnly_or_ReadsOnTargetAndJunctionCounts_set_no',
                'all_a3ss': 'A3SS_ReadsOnTargetAndJunctionCounts_event_id_set_no',
                'jcandall_a3ss': 'A3SS_JunctionCountOnly_and_ReadsOnTargetAndJunctionCounts_set_no',
                'old_a3ss': 'A3SS_old_event_id_set_no', 'jc_a3ss': 'A3SS_JunctionCountOnly_event_id_set_no',
                'novel_a3ss': 'A3SS_novel_event_id_no', 'total_a3ss': 'A3SS_all_event_id_no',
                'jcorall_a3ss': 'A3SS_JunctionCountOnly_or_ReadsOnTargetAndJunctionCounts_set_no',
                'jc_total': 'total_JunctionCountOnly_event_id_set_no', 'total_total': 'total_as_events_no',
                'all_total': 'total_ReadsOnTargetAndJunctionCounts_event_id_set_no',
                'novel_total': 'total_as_novel_events_no'}
        
        self._RMATS_PSI_MONGO_TABLE_FIELD_DIC = \
            {
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
    
    @report_check
    def add_sg_splicing_rmats(self, params=None, major=True, group=None, ref_gtf=None, name=None, outpath=None):
        task_id = self.bind_object.sheet.id
        project_sn = self.bind_object.sheet.project_sn
        if ref_gtf:
            chr_set = [e.strip() for e in
                       subprocess.check_output('awk -F \'\\t\'  \'$0!~/^#/{print $1}\' %s  | uniq | sort |uniq' % ref_gtf,
                                               shell=True).strip().split('\n')]
        else:
            raise Exception('导表时没有设置ref gtf路径')
        if not outpath:
            raise Exception('导表时没有设置rmats结果根目录路径')
        if (not group) or (not isinstance(group, dict)) or \
                (not set(group.values()) == {'s1', 's2'}) or len(group) != 2:
            raise Exception(
                '您设置的group：%s 不合法，'
                '合法的group 应该是 字典，且 '
                '固定格式为{\'a_group_name\':\'s1\', \'b_group_name\':\'s2\'}, 其中 s1和s2为固定')
        insert_data = {
            'project_sn': project_sn,
            'task_id': task_id,
            'name': name if name else 'SplicingRmats_' + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S")),
            'desc': '可变剪接rmats计算主表',
            'created_ts': datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            'params':
                json.dumps(params, sort_keys=True, separators=(',', ':')) if isinstance(params, dict) else params,
            'status': 'end',
            'group': group,
            'chr_set': chr_set,
            'rmats_out_root_dir': outpath
        }
        collection_obj = self.db['sg_splicing_rmats']
        try:
            splicing_id = collection_obj.insert_one(insert_data).inserted_id
            print('导入rmats主表成功')
        except Exception as e:
            raise Exception('导入rmats主表失败:%s' % e)
        if major:
            src_detail_file = os.path.join(outpath, 'all_events_detail_big_table.txt')
            src_stats_file = os.path.join(outpath, 'event_stats.file.txt')
            src_psi_file = os.path.join(outpath, 'psi_stats.file.txt')
            self.add_sg_splicing_rmats_detail(splicing_id=splicing_id, file=src_detail_file)
            stats_data = self._parse_dic_file(file=src_stats_file)
            psi_data = self._parse_dic_file(file=src_psi_file)
            stats_dic = self.add_sg_splicing_rmats_stats(splicing_id=splicing_id, src_dic=stats_data)
            psi_dic = self.add_sg_splicing_rmats_psi(splicing_id=splicing_id, src_dic=psi_data)
            self.add_sg_splicing_rmats_graph(splicing_id=splicing_id, stats=stats_dic, psi=psi_dic)
            # self.add_sg_splicing_rmats_model(splicing_id=splicing_id,sashimi_plot='')
    
    def add_sg_splicing_rmats_for_controller(self, splicing_id, outpath):
        src_detail_file = os.path.join(outpath, 'all_events_detail_big_table.txt')
        src_stats_file = os.path.join(outpath, 'event_stats.file.txt')
        src_psi_file = os.path.join(outpath, 'psi_stats.file.txt')
        self.add_sg_splicing_rmats_detail(splicing_id=splicing_id, file=src_detail_file)
        stats_data = self._parse_dic_file(file=src_stats_file)
        psi_data = self._parse_dic_file(file=src_psi_file)
        stats_dic = self.add_sg_splicing_rmats_stats(splicing_id=splicing_id, src_dic=stats_data)
        psi_dic = self.add_sg_splicing_rmats_psi(splicing_id=splicing_id, src_dic=psi_data)
        self.add_sg_splicing_rmats_graph(splicing_id=splicing_id, stats=stats_dic, psi=psi_dic)
    
    def _parse_dic_file(self, file):
        dic = dict(
            [(arr[0], arr[1]) for arr in [line.strip().split('\t') for line in open(file).readlines()]])
        return dic
    
    def add_sg_splicing_rmats_detail(self, splicing_id, file):
        rmats_detail_file_head_index_dic = dict(zip(self._RMATS_DETAIL_TABLE_HEAD, range(64)))
        file_tmp = os.path.join(os.path.dirname(file), os.path.basename(file) + '_tmp')
        subprocess.call('sed /^[[:space:]]*$/d %s > %s' % (file, file_tmp), shell=True)
        file_content = open(file)
        file_content.readline()
        data_lst = []
        while 1:
            data = {}
            line = file_content.readline()
            if not line.strip():
                break
            arr = re.split('\t+', line.strip())
            if len(arr) != 64:
                raise Exception(
                    'line in rmats big detail file %s is not legal: not 64 columns: %s' % (file, line))
            data['splicing_id'] = splicing_id
            for field in self._RMATS_DETAIL_MONGO_TABLE_FIELD_DIC.keys():
                data[field] = arr[rmats_detail_file_head_index_dic[self._RMATS_DETAIL_MONGO_TABLE_FIELD_DIC[field]]]

            data["pvalue_jc"] = float(data["pvalue_jc"]) if data["pvalue_jc"] != "null" else data["pvalue_jc"]
            data["pvalue_all"] = float(data["pvalue_all"]) if data["pvalue_all"] != "null" else data["pvalue_all"]
            data["fdr_jc"] = float(data["fdr_jc"]) if data["fdr_jc"] != "null" else data["fdr_jc"]
            data["fdr_all"] = float(data["fdr_all"]) if data["fdr_all"] != "null" else data["fdr_all"]
            data['splicing_id'] = splicing_id
            data['no_diff'] = 'no' if data['diff_jc_or_all'] == 'yes' else 'yes'
            data_lst.append(data)
        collection_obj = self.db['sg_splicing_rmats_detail']
        try:
            collection_obj.insert_many(data_lst)
            print("导入rmats事件详情表：%s信息成功" % (file))
        except Exception as e:
            raise Exception("导入rmats事件详情表：%s信息出错:%s" % (file, e))
        # os.r(file_tmp)
        return True
    
    def add_sg_splicing_rmats_stats(self, splicing_id, src_dic):
        tmp_dic = {}
        for field in self._RMATS_STATS_MONGO_TABLE_FIELD_DIC.keys():
            tmp_dic[field] = int(src_dic[self._RMATS_STATS_MONGO_TABLE_FIELD_DIC[field]])
        tmp_dic['old_total'] = int(tmp_dic['total_total']) - int(tmp_dic['novel_total'])
        tmp_dic['jcandall_total'] = sum(
            [tmp_dic['jcandall_mxe'], tmp_dic['jcandall_se'], tmp_dic['jcandall_ri'],
             tmp_dic['jcandall_a3ss'], tmp_dic['jcandall_a5ss']])
        tmp_dic['jcorall_total'] = sum(
            [tmp_dic['jcorall_mxe'], tmp_dic['jcorall_se'], tmp_dic['jcorall_ri'],
             tmp_dic['jcorall_a3ss'], tmp_dic['jcorall_a5ss']])
        
        as_index = dict(zip(range(3), ['old', 'novel', 'total']))
        diff_index = dict(zip(range(4), ['jc', 'all', 'jcandall', 'jcorall']))
        stats_dic = {
            'as_stats':
                {'a3ss': [0, 0, 0], 'a5ss': [0, 0, 0], 'mxe': [0, 0, 0], 'se': [0, 0, 0], 'ri': [0, 0, 0],
                 'total': [0, 0, 0]
                 },
            'diff_stats':
                {'a3ss': [0, 0, 0, 0], 'a5ss': [0, 0, 0, 0], 'mxe': [0, 0, 0, 0], 'se': [0, 0, 0, 0],
                 'ri': [0, 0, 0, 0],
                 'total': [0, 0, 0, 0]
                 },
            'splicing_id': splicing_id
        }
        
        for type in stats_dic['as_stats'].keys():
            for i in range(3):
                stats_dic['as_stats'][type][i] = int(tmp_dic[as_index[i] + '_' + type])
        for type in stats_dic['diff_stats'].keys():
            for i in range(4):
                stats_dic['diff_stats'][type][i] = int(tmp_dic[diff_index[i] + '_' + type])
        
        collection_obj = self.db['sg_splicing_rmats_stats']
        try:
            collection_obj.insert_one(stats_dic)
            print("导入rmats事件统计表信息成功")
        except Exception as e:
            raise Exception("导入rmats事件统计表信息出错 %s" % (e))
        return tmp_dic
    
    def add_sg_splicing_rmats_psi(self, splicing_id, src_dic):
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
            'splicing_id': splicing_id
        }
        psi_index = dict(zip(range(3), ['exc', 'inc', 'total']))
        for choice in ['s1_all', 's1_jc', 's2_all', 's2_jc']:
            for type in psi_dic[choice].keys():
                for i in range(3):
                    psi_dic[choice][type][i] = tmp_dic[type + '_' + choice + '_' + psi_index[i]]
        
        collection_obj = self.db['sg_splicing_rmats_psi']
        try:
            collection_obj.insert_one(psi_dic)
            print("导入rmats psi统计表信息成功")
        except Exception as e:
            raise Exception("导入rmats psi 统计表信息出错:%s" % e)
        return tmp_dic
    
    def add_sg_splicing_rmats_graph(self, splicing_id, stats, psi):
        graph_dic = {'stats_pie': {'novel': {'A3SS': 0, 'MXE': 0, 'RI': 0, 'SE': 0, 'A5SS': 0, },
                                   'old': {'A3SS': 0, 'MXE': 0, 'RI': 0, 'SE': 0, 'A5SS': 0, },
                                   'total': {'A3SS': 0, 'MXE': 0, 'RI': 0, 'SE': 0, 'A5SS': 0, }},
                     'psi_bar': {'s1': {'jc': {'exc': {'A3SS': 0, 'MXE': 0, 'RI': 0, 'SE': 0, 'A5SS': 0, },
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
                     'splicing_id': splicing_id
                     }
        
        for kind in graph_dic['stats_pie'].keys():
            for type in graph_dic['stats_pie'][kind].keys():
                graph_dic['stats_pie'][kind][type] = stats[kind + '_' + type.lower()]
        for src in graph_dic['diff_pie'].keys():
            for type in graph_dic['diff_pie'][src].keys():
                graph_dic['diff_pie'][src][type] = stats[src + '_' + type.lower()]
        for sample in graph_dic['psi_bar'].keys():
            for src in graph_dic['psi_bar'][sample].keys():
                for psi_kind in graph_dic['psi_bar'][sample][src].keys():
                    for type in graph_dic['psi_bar'][sample][src][psi_kind].keys():
                        graph_dic['psi_bar'][sample][src][psi_kind][type] = psi[
                            '_'.join([type.lower(), sample, src, psi_kind])]
        
        collection_obj = self.db['sg_splicing_rmats_graph']
        try:
            collection_obj.insert_one(graph_dic)
            print("导入rmats graph信息成功")
        except Exception as e:
            raise Exception("导入rmats graph信息出错:%s" % (e))

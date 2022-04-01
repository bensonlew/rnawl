# -*- coding: utf-8 -*-
# __author__ = fiona
# time: 2017/5/9 14:14
# last modified by shicaiping at 20180510
import re, subprocess
import json
import datetime, os
from biocluster.api.database.base import Base, report_check
from mbio.api.database.ref_rna_v2.api_base import ApiBase

class SplicingRmats(ApiBase):
    def __init__(self, bind_object):
        super(SplicingRmats, self).__init__(bind_object)
        self._project_type = 'ref_rna_v2'
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
            {
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
                'novel_total': 'total_as_novel_events_no'}

        self._RMATS_STATS_MONGO_TABLE_FIELD_DIC = \
            {
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
    def add_sg_splicing_rmats(self, params=None, major=True, group_dict=None, name=None, outpath=None, compare_plan=None):
        task_id = self.bind_object.sheet.id
        if not outpath:
            raise Exception('导表时没有设置rmats结果根目录路径')
        case_group_name = compare_plan.split('|')[0]
        control_group_name = compare_plan.split('|')[1]
        insert_data = {
            'task_id': task_id,
            'name': name if name else "SplicingRmats_" + case_group_name + "_vs_" + control_group_name + "_" + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S")),
            'desc': '可变剪接rmats计算主表',
            'created_ts': datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            'params':
                json.dumps(params, sort_keys=True, separators=(',', ':')) if isinstance(params, dict) else params,
            'status': 'end',
            'result_dir': outpath,
            'compare_plan': compare_plan
        }
        collection_obj = self.db['sg_splicing_rmats']
        try:
            splicing_id = collection_obj.insert_one(insert_data).inserted_id
            self.update_db_record('sg_splicing_rmats', splicing_id, status="end", main_id=splicing_id)
            print('导入rmats主表成功')
        except Exception as e:
            raise Exception('导入rmats主表失败:%s' % e)
        if major:
            src_detail_file = os.path.join(outpath, 'all_events_detail_big_table.txt')
            type_stats_file = os.path.join(outpath, 'event_type.file.txt')

            self.add_sg_splicing_rmats_detail(splicing_id=splicing_id, file=src_detail_file)
            type_stats_data = self._parse_dic_file(file=type_stats_file)
            self.add_sg_splicing_rmats_type_stats(splicing_id=splicing_id, src_dic=type_stats_data)
            self.add_sg_splicing_rmats_stats(splicing_id=splicing_id, major=True, outpath=outpath)
    
    def add_sg_splicing_rmats_for_controller(self, splicing_id, outpath):
        src_detail_file = os.path.join(outpath, 'all_events_detail_big_table.txt')
        type_stats_file = os.path.join(outpath, 'event_type.file.txt')

        self.add_sg_splicing_rmats_detail(splicing_id=splicing_id, file=src_detail_file)
        type_stats_data = self._parse_dic_file(file=type_stats_file)
        self.add_sg_splicing_rmats_type_stats(splicing_id=splicing_id, src_dic=type_stats_data)
        self.add_sg_splicing_rmats_stats(splicing_id=splicing_id, major=True, outpath=outpath)

    def add_sg_rmats_stat_for_controller(self, stat_id, outpath):
        diff_stats_file = os.path.join(outpath, 'event_stats.file.txt')
        psi_stat_file = os.path.join(outpath, 'psi_stats.file.txt')

        psi_data = self._parse_dic_file(file=psi_stat_file)
        diff_stats_data = self._parse_dic_file(file=diff_stats_file)
        diff_stats_dic= self.add_sg_splicing_rmats_diff_stats(stat_id=stat_id, src_dic=diff_stats_data)
        psi_dic = self.add_sg_splicing_rmats_psi(stat_id=stat_id, src_dic=psi_data)
        self.add_sg_splicing_rmats_graph(stat_id=stat_id, stats=diff_stats_dic, psi=psi_dic)
    
    def _parse_dic_file(self, file):
        dic = dict(
            [(arr[0], arr[1]) for arr in [line.strip().split('\t') for line in open(file).readlines()]])
        return dic
    
    def add_sg_splicing_rmats_detail(self, splicing_id, file):
        rmats_detail_file_head_index_dic = dict(zip(self._RMATS_DETAIL_TABLE_HEAD, range(64)))
        #file_tmp = os.path.join(os.path.dirname(file), os.path.basename(file) + '_tmp')
        subprocess.call('sed -i /^[[:space:]]*$/d %s' % (file), shell=True)
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
            data_lst.append(data)
        collection_obj = self.db['sg_splicing_rmats_detail']
        file_content = open(file)
        if len(file_content.readlines()) >= 2:
            try:
                collection_obj.insert_many(data_lst)
                print("导入rmats事件详情表：%s信息成功" % (file))
            except Exception as e:
                raise Exception("导入rmats事件详情表：%s信息出错:%s" % (file, e))
            # os.r(file_tmp)
        return True
    
    def add_sg_splicing_rmats_stats(self, params=None, major=True, splicing_id=None, outpath=None, name=None):
        if not outpath:
            raise Exception('导表时没有设置差异可变剪接事件统计结果根目录路径')
        if major:
            task_id = self.db['sg_splicing_rmats'].find_one({"_id": splicing_id})["task_id"]
            compare_plan = self.db['sg_splicing_rmats'].find_one({"_id": splicing_id})['compare_plan']
            case_group_name = compare_plan.split('|')[0]
            control_group_name = compare_plan.split('|')[1]
            group = {case_group_name: 's1', control_group_name: 's2'}
            params = {'splicing_id': str(splicing_id), 'pvalue_fdr': "fdr", 'submit_location' : 'splicingrmats_type',
                    'fdr': "0.05", 'task_id': task_id, 'task_type': "2", 'psi': "0"}
            insert_data = {
                'task_id': task_id,
                'name': name if name else "RmatsStat_" + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S")),
                'desc': '差异可变剪接事件统计主表',
                'created_ts': datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
                'params':
                    json.dumps(params, sort_keys=True, separators=(',', ':')) if isinstance(params, dict) else params,
                'status': 'end',
                'group': group
            }
            collection_obj = self.db['sg_splicing_rmats_stats']
            try:
                stat_id = collection_obj.insert_one(insert_data).inserted_id
                self.update_db_record('sg_splicing_rmats_stats', stat_id, status="end", main_id=stat_id)
                print('导入差异可变剪接事件统计主表成功')
            except Exception as e:
                raise Exception('导入差异可变剪接事件统计主表失败:%s' % e)

            diff_stats_file = os.path.join(outpath, 'event_stats.file.txt')
            psi_stat_file = os.path.join(outpath, 'psi_stats.file.txt')

            psi_data = self._parse_dic_file(file=psi_stat_file)
            diff_stats_data = self._parse_dic_file(file=diff_stats_file)
            diff_stats_dic = self.add_sg_splicing_rmats_diff_stats(stat_id=stat_id, src_dic=diff_stats_data)
            psi_dic = self.add_sg_splicing_rmats_psi(stat_id=stat_id, src_dic=psi_data)
            self.add_sg_splicing_rmats_graph(stat_id=stat_id, stats=diff_stats_dic, psi=psi_dic)

    def add_sg_splicing_rmats_type_stats(self, splicing_id, src_dic):
        tmp_dic = {}
        for field in self._RMATS_TYPE_MONGO_TABLE_FIELD_DIC.keys():
            tmp_dic[field] = int(src_dic[self._RMATS_TYPE_MONGO_TABLE_FIELD_DIC[field]])
        tmp_dic['old_total'] = int(tmp_dic['total_total']) - int(tmp_dic['novel_total'])
        as_index = dict(zip(range(3), ['old', 'novel', 'total']))

        stats_dic = {
            'as_stats':
                {'a3ss': [0, 0, 0], 'a5ss': [0, 0, 0], 'mxe': [0, 0, 0], 'se': [0, 0, 0], 'ri': [0, 0, 0],
                 'total': [0, 0, 0]
                 },
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
        try:
            collection_obj.insert_one(stats_dic)
            print("导入rmats事件类型统计表信息成功")
        except Exception as e:
            raise Exception("导入rmats事件类型统计表信息出错 %s" % (e))
        return tmp_dic

    def add_sg_splicing_rmats_diff_stats(self, stat_id, src_dic):
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
            'diff_stats':
                {'a3ss': [0, 0, 0, 0], 'a5ss': [0, 0, 0, 0], 'mxe': [0, 0, 0, 0], 'se': [0, 0, 0, 0],
                 'ri': [0, 0, 0, 0],
                 'total': [0, 0, 0, 0]
                 },
            'stat_id': stat_id,
        }

        for type in stats_dic['diff_stats'].keys():
            for i in range(4):
                stats_dic['diff_stats'][type][i] = int(tmp_dic[diff_index[i] + '_' + type])

        collection_obj = self.db['sg_splicing_rmats_diff_stats']
        try:
            main_id = collection_obj.insert_one(stats_dic).inserted_id
            print("导入rmats事件统计表信息成功")
        except Exception as e:
            raise Exception("导入rmats事件统计表信息出错 %s" % (e))
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
        try:
            collection_obj.insert_one(psi_dic)
            print("导入rmats psi统计表信息成功")
        except Exception as e:
            raise Exception("导入rmats psi 统计表信息出错:%s" % e)
        return tmp_dic
    
    def add_sg_splicing_rmats_graph(self, stat_id, stats, psi):
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
                     'stat_id': stat_id
                     }

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

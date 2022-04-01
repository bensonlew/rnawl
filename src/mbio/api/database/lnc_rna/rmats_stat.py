# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

from mbio.api.database.lnc_rna.api_base import ApiBase
from biocluster.api.database.base import report_check
import os
from bson.objectid import ObjectId
import unittest

class RmatsStat(ApiBase):
    '''
    last_modify: 2019.04.02
    '''
    def __init__(self, bind_object):
        super(RmatsStat, self).__init__(bind_object)
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

    def parse_dic_file(self, file):
        return dict([(arr[0], arr[1]) for arr in [line.strip().split('\t') for line in open(file).readlines()]])

    @report_check
    def add_rmats_stat_detail(self, outpath, main_id):
        main_id = ObjectId(main_id)
        psi_data = self.parse_dic_file(file=os.path.join(outpath, 'psi_stats.file.txt'))
        event_data = self.parse_dic_file(file=os.path.join(outpath, 'event_stats.file.txt'))
        diff_stats_dic = self.add_sg_splicing_rmats_diff_stats(stat_id=main_id, src_dic=event_data, outpath=outpath)
        psi_dic = self.add_sg_splicing_rmats_psi(stat_id=main_id, src_dic=psi_data)
        ret = self.add_sg_splicing_rmats_graph(stat_id=main_id, stats=diff_stats_dic, psi=psi_dic, outpath=outpath)
        if ret:
            insert_dict = {
                'main_id': ObjectId(main_id),
                'status': 'end'
            }
            self.update_db_record('sg_splicing_rmats_stats', main_id, insert_dict=insert_dict)

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

        # export to excel on 20201215
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

class TestFunction(unittest.TestCase):
    '''
    This is test for the api. Just run this script to do test.
    '''
    def test(self):
        import random
        from mbio.workflows.lnc_rna.lnc_rna_test_api import LncRnaTestApiWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'rmats_model_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
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
        wf.test_api = wf.api.api('lnc_rna.rmats_stat')
        wf.test_api.add_rmats_stat_detail(
            outpath='/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/lnc_rna/rmats_stat/output',
            main_id='5c92247417b2bf779611497e'
        )

if __name__ == '__main__':
    unittest.main()
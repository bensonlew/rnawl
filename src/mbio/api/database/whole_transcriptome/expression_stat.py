# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import datetime
import json
import unittest

import pandas as pd

from mbio.api.database.whole_transcriptome.api_base import ApiBase


class ExpressionStat(ApiBase):
    def __init__(self, bind_object):
        super(ExpressionStat, self).__init__(bind_object)

    def add_expression_stat(self, map_dict, task_id, project_sn):
        time_now = datetime.datetime.now()
        name = 'expression_stat_{}'.format(time_now.strftime('%Y%m%d_%H%M%S'))
        params = json.dumps({'task_id': task_id, 'submit_location': 'expressionstat', 'task_type': 2}, sort_keys=True)
        main_dict = {
            'task_id': task_id,
            'project_sn': project_sn,
            'name': name,
            'created_ts': time_now.strftime('%Y-%m-%d %H:%M:%S'),
            'desc': 'Workflow result',
            'params': params,
            'status': 'start',
            'version': 'v1'
        }
        main_id = self.create_db_table('expression_stat', [main_dict])
        c_count_matrix = map_dict['c_count'] if 'c_count' in map_dict else None
        self.add_expression_stat_detail_t(main_id, map_dict['t_type'], map_dict['t_count'], map_dict['rna_select'],
                                          c_count_matrix)
        self.add_expression_stat_detail_g(main_id, map_dict['g_type'], map_dict['g_count'], map_dict['rna_select'])
        self.update_db_record('expression_stat', main_id, insert_dict={'main_id': main_id, 'status': 'end'})

    def add_expression_stat_detail_t(self, stat_id, type_table, t_count_matrix, t_select, c_count_matrix=None):
        self.bind_object.logger.info('start inserting documents into mongo with ({})'.format({'level': 'T'}))
        m_dict = {'category': 'mRNA', 'known': 0, 'novel': 0, 'total': 0}
        l_dict = {'category': 'lncRNA', 'known': 0, 'novel': 0, 'total': 0}
        type_df = pd.read_table(type_table, names=['transcript_id', 'gene_id', 'category', 'kind'], index_col=0)
        longrna_ids = list(pd.read_table(t_count_matrix)['seq_id'].unique())
        has_lnc = True if 'lncRNA' in t_select else False
        for longrna_id in longrna_ids:
            type_series = type_df.loc[longrna_id]
            if type_series['category'] == 'mRNA':
                m_dict[type_series['kind']] += 1
                m_dict['total'] += 1
            elif type_series['category'] == 'lncRNA' and has_lnc:
                l_dict[type_series['kind']] += 1
                l_dict['total'] += 1
        else:
            data = [m_dict, l_dict]
        if c_count_matrix:
            c_dict = {'category': 'circRNA', 'known': 0,
                      'novel': len(pd.read_table(c_count_matrix)['seq_id'].unique()),
                      'total': len(pd.read_table(c_count_matrix)['seq_id'].unique())}
            data.append(c_dict)
        self.create_db_table('expression_stat_detail', data, {'level': 'T', 'stat_id': stat_id})

    def add_expression_stat_detail_g(self, stat_id, type_table, g_count_matrix, g_select):
        self.bind_object.logger.info('start inserting documents into mongo with ({})'.format({'level': 'G'}))
        m_dict = {'category': 'mRNA', 'known': 0, 'novel': 0, 'total': 0}
        l_dict = {'category': 'lncRNA', 'known': 0, 'novel': 0, 'total': 0}
        type_df = pd.read_table(type_table, names=['gene_id', 'transcript_ids', 'category', 'kind'], index_col=0)
        longrna_gene_ids = list(pd.read_table(g_count_matrix)['seq_id'].unique())
        has_lnc = True if 'lncRNA' in g_select else False
        for gene_id in longrna_gene_ids:
            type_series = type_df.loc[gene_id]
            if type_series['category'] == 'mRNA':
                m_dict[type_series['kind']] += 1
                m_dict['total'] += 1
            elif type_series['category'] == 'lncRNA' and has_lnc:
                l_dict[type_series['kind']] += 1
                l_dict['total'] += 1
        else:
            data = [m_dict, l_dict]
        self.create_db_table('expression_stat_detail', data, {'level': 'G', 'stat_id': stat_id})


class TestFunction(unittest.TestCase):
    '''
    This is test for the api. Just run this script to do test.
    '''

    def test(self):
        import random
        from mbio.workflows.whole_transcriptome.whole_transcriptome_test_api import WholeTranscriptomeTestApiWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'expression_stat_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'workflow',
            'name': 'whole_transcriptome.whole_transcriptome_test_api',
            'options': {}
        }
        wheet = Sheet(data=data)
        wf = WholeTranscriptomeTestApiWorkflow(wheet)
        wf.IMPORT_REPORT_DATA = True
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.test_api = wf.api.api('whole_transcriptome.expression_stat')
        map_dict = {
            't_type': '/mnt/ilustre/users/sanger-dev/workspace/20191014/Longrna_workflow_6700_7886/LargeGush/output/filter_by_express/filtered_file/trans_type.xls',
            't_count': '/mnt/ilustre/users/sanger-dev/workspace/20191015/Single_exp_make_6370_6765/ExpMake/output/count/T.reads.txt',
            'c_count': '/mnt/ilustre/users/sanger-dev/workspace/20191015/Single_exp_make_6370_6765/ExpMake/output/count/C.reads.txt',
            'g_type': '/mnt/ilustre/users/sanger-dev/workspace/20191014/Longrna_workflow_6700_7886/LargeGush/output/filter_by_express/filtered_file/gene_type.xls',
            'g_count': '/mnt/ilustre/users/sanger-dev/workspace/20191015/Single_exp_make_6370_6765/ExpMake/output/count/G.reads.txt'
        }
        wf.test_api.add_expression_stat(
            map_dict=map_dict,
            task_id='whole_transcriptome',
            project_sn='whole_transcriptome'
        )


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)

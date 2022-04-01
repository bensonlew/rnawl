# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import datetime
import json
import os
import unittest

import pandas as pd

from mbio.api.database.whole_transcriptome.api_base import ApiBase


class DiffExpStat(ApiBase):
    def __init__(self, bind_object):
        super(DiffExpStat, self).__init__(bind_object)

    def add_diff_exp_stat(self, map_dict, arg_dict, task_id, project_sn, library):
        time_now = datetime.datetime.now()
        name = 'diff_exp_stat_{}'.format(time_now.strftime('%Y%m%d_%H%M%S'))
        params = json.dumps({'task_id': task_id, 'submit_location': 'diffexpstat', 'task_type': 2}, sort_keys=True)
        main_dict = {
            'task_id': task_id,
            'project_sn': project_sn,
            'name': name,
            'created_ts': time_now.strftime('%Y-%m-%d %H:%M:%S'),
            'desc': 'Workflow result',
            'params': params,
            'status': 'start',
            'library': library,
            'version': 'v1'
        }
        for k, v in arg_dict.items():
            main_dict.update({k: v})
        main_id = self.create_db_table('diff_exp_stat', [main_dict])
        if library == 'long':
            c_output_dir = map_dict['c_output_dir'] if 'c_output_dir' in map_dict else None
            self.add_diff_exp_stat_detail_l(main_id, map_dict['control'], map_dict['t_type'], map_dict['t_output_dir'],
                                            c_output_dir)
        elif library == 'small':
            self.add_diff_exp_stat_detail_s(main_id, map_dict['control'], map_dict['s_output_dir'])
        self.update_db_record('diff_exp_stat', main_id, insert_dict={'main_id': main_id, 'status': 'end'})

    def add_diff_exp_stat_detail_l(self, stat_id, control_file, type_table, t_output_dir, c_output_dir=None):
        self.bind_object.logger.info('start inserting documents into mongo with ({})'.format({'library': 'long'}))
        cmp_dict = dict()
        for line in open(control_file):
            if line.strip() and line[0] != '#':
                ctrl, case = line.strip().split('\t')[:2]
                cmp_dict['{}_vs_{}'.format(ctrl, case)] = {rna: [0] * 3 for rna in ('mrna', 'lncrna', 'circrna')}
        type_df = pd.read_table(type_table, names=['transcript_id', 'gene_id', 'category', 'kind'], index_col=0)
        for compare in cmp_dict:
            t_detail_df = pd.read_table(os.path.join(t_output_dir, '{}.detail.txt'.format(compare)))
            for i, row in t_detail_df.iterrows():
                seq_id = row['seq_id']
                category = type_df.loc[seq_id, 'category'].lower()
                if row['significant'] == 'yes':
                    cmp_dict[compare][category][0] += 1
                    if row['regulate'] == 'up':
                        cmp_dict[compare][category][1] += 1
                    elif row['regulate'] == 'down':
                        cmp_dict[compare][category][2] += 1
            if c_output_dir:
                c_detail_df = pd.read_table(os.path.join(c_output_dir, '{}.detail.txt'.format(compare)))
                for i, row in c_detail_df.iterrows():
                    seq_id = row['seq_id']
                    if row['significant'] == 'yes':
                        cmp_dict[compare]['circrna'][0] += 1
                        if row['regulate'] == 'up':
                            cmp_dict[compare]['circrna'][1] += 1
                        elif row['regulate'] == 'down':
                            cmp_dict[compare]['circrna'][2] += 1
        data = list()
        for compare in cmp_dict:
            document = {'diff_group': compare}
            for rna, sig_nums in cmp_dict[compare].items():
                document[rna] = '|'.join(map(str, sig_nums))
            else:
                data.append(document)
        else:
            self.create_db_table('diff_exp_stat_detail', data, {'level': 'T', 'stat_id': stat_id})

    def add_diff_exp_stat_detail_s(self, stat_id, control_file, s_output_dir):
        self.bind_object.logger.info('start inserting documents into mongo with ({})'.format({'library': 'small'}))
        cmp_dict = dict()
        for line in open(control_file):
            if line.strip() and line[0] != '#':
                ctrl, case = line.strip().split('\t')[:2]
                cmp_dict['{}_vs_{}'.format(ctrl, case)] = {'up': 0, 'down': 0, 'total': 0}
        for compare in cmp_dict:
            detail_df = pd.read_table(os.path.join(s_output_dir, '{}.detail.txt'.format(compare)))
            for i, row in detail_df.iterrows():
                if row['significant'] == 'yes':
                    cmp_dict[compare]['total'] += 1
                    cmp_dict[compare][row['regulate']] += 1
        data = list()
        for compare, value in cmp_dict.items():
            document = {'diff_group': compare}
            document.update(value)
            data.append(document)
        else:
            self.create_db_table('diff_exp_stat_detail', data, {'stat_id': stat_id})


class TestFunction(unittest.TestCase):
    '''
    This is test for the api. Just run this script to do test.
    '''

    def test(self):
        import random
        from mbio.workflows.whole_transcriptome.whole_transcriptome_test_api import WholeTranscriptomeTestApiWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'diff_exp_stat_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'workflow',
            'name': 'whole_transcriptome.whole_transcriptome_test_api',
            'options': {}
        }
        wheet = Sheet(data=data)
        wf = WholeTranscriptomeTestApiWorkflow(wheet)
        wf.IMPORT_REPORT_DATA = True
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.test_api = wf.api.api('whole_transcriptome.diff_exp_stat')
        map_dict = {
            'control': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/whole_transcriptome/MJ20190118060/largeRNA/example_control.txt',
            't_type': '/mnt/ilustre/users/sanger-dev/workspace/20191016/Longrna_workflow_8323_6543/LargeGush/output/filter_by_express/filtered_file/trans_type.xls',
            't_output_dir': '/mnt/ilustre/users/sanger-dev/workspace/20191017/Single_diff_exp_5116_1524/DiffExp/output',
            'c_output_dir': '/mnt/ilustre/users/sanger-dev/workspace/20191017/Single_diff_exp_1583_7158/DiffExp/output',
        }
        wf.test_api.add_diff_exp_stat(
            map_dict=map_dict,
            task_id='whole_transcriptome',
            project_sn='whole_transcriptome'
        )


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)

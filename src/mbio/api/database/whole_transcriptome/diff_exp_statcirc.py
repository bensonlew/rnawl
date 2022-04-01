# -*- coding: utf-8 -*-
# __author__ = 'zoujiaxun'

import datetime
import json
import os
import unittest

import pandas as pd

from mbio.api.database.whole_transcriptome.api_base import ApiBase


class DiffExpStatcirc(ApiBase):
    def __init__(self, bind_object):
        super(DiffExpStatcirc, self).__init__(bind_object)

    def add_diff_exp_statcirc(self, map_dict, arg_dict, task_id, project_sn, library):
        time_now = datetime.datetime.now()
        name = 'diff_exp_statcirc_{}'.format(time_now.strftime('%Y%m%d_%H%M%S'))
        params = json.dumps({'task_id': task_id, 'submit_location': 'diffexpstat', 'task_type': 2}, sort_keys=True)
        main_dict = {
            'task_id': task_id,
            'project_sn': project_sn,
            'name': name,
            'created_ts': time_now.strftime('%Y-%m-%d %H:%M:%S'),
            'desc': 'Workflow result',
            'params': params,
            'status': 'start',
            'library': library

        }
        for k, v in arg_dict.items():
            main_dict.update({k: v})
        main_id = self.create_db_table('diff_exp_stat', [main_dict])
        # if library == 'long':
        #     c_output_dir = map_dict['c_output_dir'] if 'c_output_dir' in map_dict else None
        #     self.add_diff_exp_stat_detail_l(main_id, map_dict['control'], map_dict['t_type'], map_dict['t_output_dir'],
        #                                     c_output_dir)
        # elif library == 'small':
        self.add_diff_exp_stat_detail_c(main_id, map_dict['control'], map_dict['c_output_dir'])
        self.update_db_record('diff_exp_stat', main_id, insert_dict={'main_id': main_id, 'status': 'end'})

    def add_diff_exp_stat_detail_c(self, stat_id, control_file, c_output_dir):
        self.bind_object.logger.info('start inserting documents into mongo with ({})'.format({'library': 'circ'}))
        cmp_dict = dict()
        for line in open(control_file):
            if line.strip() and line[0] != '#':
                ctrl, case = line.strip().split('\t')[:2]
                cmp_dict['{}_vs_{}'.format(ctrl, case)] = {'up': 0, 'down': 0, 'total': 0}
        for compare in cmp_dict:
            detail_df = pd.read_table(os.path.join(c_output_dir, '{}.detail.txt'.format(compare)))
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
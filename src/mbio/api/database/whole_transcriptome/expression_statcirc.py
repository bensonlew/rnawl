# -*- coding: utf-8 -*-
# __author__ = 'zoujiaxun'

import datetime
import json
import unittest

import pandas as pd

from mbio.api.database.whole_transcriptome.api_base import ApiBase


class ExpressionStatcirc(ApiBase):
    def __init__(self, bind_object):
        super(ExpressionStatcirc, self).__init__(bind_object)

    def add_expression_statcirc(self, stat, task_id, project_sn, group=None):
        if group:
            sample_list = list()
            with open(group, 'r') as g:
                for line in g.readlines():
                    if line.startswith('#'):
                        continue
                    sample_list.append(line.strip().split('\t')[0])
            sample_list.append('total')
        time_now = datetime.datetime.now()
        name = 'expression_stat_{}'.format(time_now.strftime('%Y%m%d_%H%M%S'))
        params = json.dumps({'task_id': task_id, 'submit_location': 'expressionstat', 'task_type': 2}, sort_keys=True)
        df = pd.read_table(stat)
        main_dict = {
            'task_id': task_id,
            'project_sn': project_sn,
            'name': name,
            'created_ts': time_now.strftime('%Y-%m-%d %H:%M:%S'),
            'desc': 'Workflow result',
            'params': params,
            'status': 'start'
        }
        main_id = self.create_db_table('expression_stat', [main_dict])
        df['stat_id'] = main_id
        data_list = df.to_dict('r')
        if group:
            data_list.sort(key=lambda x: sample_list.index(x['sample']))
        self.create_db_table('expression_stat_detail', data_list)
        self.update_db_record('expression_stat', main_id, insert_dict={'main_id': main_id, 'status': 'end'})

class TestFunction(unittest.TestCase):
    '''
    This is test for the api. Just run this script to do test.
    '''

    def test(self):
        import random
        from mbio.workflows.whole_transcriptome.whole_transcriptome_test_api import WholeTranscriptomeTestApiWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'exp_stat_circ_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'workflow',
            'name': 'whole_transcriptome.whole_transcriptome_test_api',
            'options': {}
        }
        wheet = Sheet(data=data)
        wf = WholeTranscriptomeTestApiWorkflow(wheet)
        wf.IMPORT_REPORT_DATA = True
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.test_api = wf.api.api('whole_transcriptome.expression_statcirc')

        wf.test_api.add_expression_statcirc(
            task_id='whole_transcriptome',
            project_sn='whole_transcriptome',
            stat= '/mnt/ilustre/users/sanger-dev/workspace/20191125/Single_exp_make_5073_4969/ExpMakecirc/output/count/C.stat.txt'
        )


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
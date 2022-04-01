# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import datetime
import json
import os
import unittest

import pandas as pd
from biocluster.api.database.base import report_check

from mbio.api.database.whole_transcriptome.api_base import ApiBase


class RmatsCount(ApiBase):
    def __init__(self, bind_object):
        super(RmatsCount, self).__init__(bind_object)

    @report_check
    def add_rmats_count(self, outpath, group=None):
        if group:
            sample_list = list()
            with open(group, 'r') as g:
                for line in g.readlines():
                    if line.startswith('#'):
                        continue
                    sample_list.append(line.strip().split('\t')[0])
        jc_df = pd.read_table(os.path.join(outpath, 'sample.event.count.JC.txt'))
        jc_df.sort_values('SAMPLE', inplace=True)
        jc_df = jc_df.rename(columns=str.lower)
        jcec_df = pd.read_table(os.path.join(outpath, 'sample.event.count.JCEC.txt'))
        jcec_df.sort_values('SAMPLE', inplace=True)
        jcec_df = jcec_df.rename(columns=str.lower)
        task_id = self.bind_object.sheet.id
        project_sn = self.bind_object.sheet.project_sn
        time_now = datetime.datetime.now()
        name = 'Splicing_count_{}'.format(time_now.strftime('%Y%m%d_%H%M%S'))
        created_ts = time_now.strftime('%Y-%m-%d %H:%M:%S')
        params = json.dumps({
            'task_id': task_id,
            'submit_location': 'splicingrmats_count',
            'task_type': 2
        })
        main_info = {
            'task_id': task_id,
            'project_sn': project_sn,
            'name': name,
            'desc': 'alternative splicing count main table',
            'created_ts': time_now.strftime('%Y-%m-%d %H:%M:%S'),
            'params': params,
            'status': 'start'
        }
        main_id = self.create_db_table('splicing_rmats_count', [main_info])
        data_list_jc = jc_df.to_dict('r')
        data_list_jcec = jcec_df.to_dict('r')
        if group:
            data_list_jc.sort(key=lambda x: sample_list.index(x['sample']))
            data_list_jcec.sort(key=lambda x: sample_list.index(x['sample']))
        # self.create_db_table('splicing_rmats_count_detail', jc_df.to_dict('r'),
        #                      {'as_type': 'JC', 'count_id': main_id})
        # self.create_db_table('splicing_rmats_count_detail', jcec_df.to_dict('r'),
        #                      {'as_type': 'JCEC', 'count_id': main_id})
        self.create_db_table('splicing_rmats_count_detail', data_list_jc,
                             {'as_type': 'JC', 'count_id': main_id})
        self.create_db_table('splicing_rmats_count_detail', data_list_jcec,
                             {'as_type': 'JCEC', 'count_id': main_id})
        self.update_db_record('splicing_rmats_count', main_id, insert_dict={'main_id': main_id, 'status': 'end'})


class TestFunction(unittest.TestCase):
    '''
    This is test for the api. Just run this script to do test.
    '''

    def test(self):
        import random
        from mbio.workflows.whole_transcriptome.whole_transcriptome_test_api import WholeTranscriptomeTestApiWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'rmats_count_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'workflow',
            'name': 'whole_transcriptome.whole_transcriptome_test_api',
            'options': {}
        }
        wheet = Sheet(data=data)
        wf = WholeTranscriptomeTestApiWorkflow(wheet)
        wf.sheet.id = 'whole_transcriptome'
        wf.sheet.project_sn = 'whole_transcriptome'
        wf.IMPORT_REPORT_DATA = True
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.test_api = wf.api.api('whole_transcriptome.rmats_count')
        wf.test_api.add_rmats_count(
            outpath='/mnt/ilustre/users/sanger-dev/workspace/20190717/Refrna_tsg_34883/Rmats/output'
        )


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)

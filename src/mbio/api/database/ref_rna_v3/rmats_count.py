# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

from mbio.api.database.ref_rna_v2.api_base import ApiBase
from biocluster.api.database.base import report_check
import pandas as pd
import datetime
import json
import os
from bson.objectid import ObjectId
import unittest

class RmatsCount(ApiBase):
    def __init__(self, bind_object):
        super(RmatsCount, self).__init__(bind_object)

    @report_check
    def add_rmats_count(self, outpath, group):
        sample_list = list()
        with open(group, 'r') as g:
            for line in g.readlines():
                if line.startswith('#'):
                    continue
                sample_list.append(line.strip().split('\t')[0])
        jc_df = pd.read_table(os.path.join(outpath, 'sample.event.count.JC.txt'))
        jc_df = jc_df.rename(columns=str.lower)
        samples = list(jc_df["sample"])
        jcec_df = pd.read_table(os.path.join(outpath, 'sample.event.count.JCEC.txt'))
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
        main_id = self.create_db_table('sg_splicing_rmats_count', [main_info])
        jc_df_list = jc_df.to_dict('r')
        jc_df_list.sort(key=lambda x: sample_list.index(x['sample']))
        jcec_df_list = jcec_df.to_dict('r')
        jcec_df_list.sort(key=lambda x: sample_list.index(x['sample']))
        self.create_db_table('sg_splicing_rmats_count_detail', jc_df_list,
                             {'as_type': 'JC', 'count_id': main_id})
        self.create_db_table('sg_splicing_rmats_count_detail', jcec_df_list,
                             {'as_type': 'JCEC', 'count_id': main_id})
        self.update_db_record('sg_splicing_rmats_count', main_id, insert_dict={'main_id': main_id, 'status': 'end', 'samples': samples})

class TestFunction(unittest.TestCase):
    '''
    This is test for the api. Just run this script to do test.
    '''
    def test(self):
        import random
        from mbio.workflows.ref_rna_v2.refrna_test_api import RefrnaTestApiWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'rmats_count_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'workflow',
            'name': 'ref_rna_v2.refrna_test_api',
            'options': {}
        }
        wheet = Sheet(data=data)
        wf = RefrnaTestApiWorkflow(wheet)
        wf.sheet.id = 'ref_rna_v3'
        wf.sheet.project_sn = 'ref_rna_v3'
        wf.IMPORT_REPORT_DATA = True
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.test_api = wf.api.api('ref_rna_v3.rmats_count')
        wf.test_api.add_rmats_count(
            outpath='/mnt/ilustre/users/sanger-dev/workspace/20190618/Single_rmats_6811_9583/Rmats/output'
        )

if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)

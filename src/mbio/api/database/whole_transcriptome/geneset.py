# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import datetime
import json
import os
import unittest

import pandas as pd

from mbio.api.database.whole_transcriptome.api_base import ApiBase


class Geneset(ApiBase):
    def __init__(self, bind_object):
        super(Geneset, self).__init__(bind_object)

    def add_geneset_diff(self, map_dict, level, category, source, task_id, project_sn):
        for vs_pair, detail_table in map_dict.items():
            df = pd.read_table(detail_table)
            if 'category' not in df:
                df['category'] = category
            df = df[df['significant'] == 'yes']
            length = len(df)
            if not length:
                continue
            time_now = datetime.datetime.now()
            name = '{}_{}_{}'.format(level, category, vs_pair)
            params = json.dumps({'task_id': task_id, 'submit_location': 'genesetdiff', 'task_type': 2},
                                sort_keys=True, separators=(',', ':'))
            main_dict = {'task_id': task_id,
                         'project_sn': project_sn,
                         'name': name,
                         'created_ts': time_now.strftime('%Y-%m-%d %H:%M:%S'),
                         'params': params,
                         'status': 'start',
                         'level': level,
                         'length': length,
                         'source': source,
                         'desc': str(),
                         'type': '{}:{}'.format(category, length),
                         'is_use': 0,
                         'version': 'v1'}
            main_id = self.create_db_table('geneset', [main_dict])
            arg_dict = {'seq_list': df['seq_id'].tolist(),
                        'category_list': df['category'].tolist(),
                        'kind_list': df['kind'].tolist(),
                        'regulate_list': df['regulate'].tolist()}
            self.add_geneset_detail(arg_dict, main_id)
            self.update_db_record('geneset', main_id, insert_dict={'main_id': main_id, 'status': 'end'})

    def add_geneset_detail(self, arg_dict, geneset_id):
        self.create_db_table('geneset_detail', [arg_dict], {'geneset_id': geneset_id})


class TestFunction(unittest.TestCase):
    """
    This is test for the api. Just run this script to do test.
    """

    def test(self):
        import random
        from mbio.workflows.whole_transcriptome.whole_transcriptome_test_api import WholeTranscriptomeTestApiWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'geneset_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'workflow',
            'name': 'whole_transcriptome.whole_transcriptome_test_api',
            'options': {}
        }
        wheet = Sheet(data=data)
        wf = WholeTranscriptomeTestApiWorkflow(wheet)
        wf.IMPORT_REPORT_DATA = True
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.test_api = wf.api.api('whole_transcriptome.geneset')
        map_dict = dict()
        diff_dir = '/mnt/ilustre/users/sanger-dev/workspace/20191106/WholeTranscriptome_tsg_36088/DiffSplit/output/mrna'
        for fname in os.listdir(diff_dir):
            if fname.endswith('.detail.txt'):
                map_dict[fname[:-11]] = os.path.join(diff_dir, fname)
        wf.test_api.add_geneset_diff(
            map_dict=map_dict,
            level='T',
            category='mRNA',
            source='DE_mR_detail',
            task_id='tsg_36088',
            project_sn='188_5dba6f542345b'
        )


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)

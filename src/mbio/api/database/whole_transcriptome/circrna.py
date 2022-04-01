# -*- coding: utf-8 -*-
# __author__ = 'zoujiaxun'

import datetime
import json
import unittest

import pandas as pd

from mbio.api.database.whole_transcriptome.api_base import ApiBase


class Circrna(ApiBase):
    def __init__(self, bind_object):
        super(Circrna, self).__init__(bind_object)

    def add_circrna(self, task_id, project_sn, circ_detail, taxonomy):
        time_now = datetime.datetime.now()
        name = 'CircStat_{}'.format(time_now.strftime('%Y%m%d_%H%M%S'))
        params = json.dumps({'task_id': task_id, 'submit_location': 'circrna', 'task_type': 2}, sort_keys=True)
        df = pd.read_table(circ_detail)
        df_columns = list(df)
        if 'circbase' in df_columns:
            if taxonomy == 'Animal':
                main_dict = {
                    'task_id': task_id,
                    'project_sn': project_sn,
                    'name': name,
                    'created_ts': time_now.strftime('%Y-%m-%d %H:%M:%S'),
                    'desc': 'Circrna prediction main table',
                    'params': params,
                    'status': 'start',
                    'circbase': True,
                    'circbase_url': 'http://www.circbase.org/',
                    'type': dict(df['circrna_type'].value_counts())
                }
            if taxonomy == 'Plant':
                main_dict = {
                    'task_id': task_id,
                    'project_sn': project_sn,
                    'name': name,
                    'created_ts': time_now.strftime('%Y-%m-%d %H:%M:%S'),
                    'desc': 'Circrna prediction main table',
                    'params': params,
                    'status': 'start',
                    'circbase': True,
                    'circbase_url': 'http://ibi.zju.edu.cn/plantcircbase/',
                    'type': dict(df['circrna_type'].value_counts())
                }
        else:
            main_dict = {
                'task_id': task_id,
                'project_sn': project_sn,
                'name': name,
                'created_ts': time_now.strftime('%Y-%m-%d %H:%M:%S'),
                'desc': 'Circrna prediction main table',
                'params': params,
                'status': 'start',
                'circbase': False,
                'type': dict(df['circrna_type'].value_counts())
            }
        main_dict['version'] = 'v1'
        main_id = self.create_db_table('circrna_identify', [main_dict])
        df['detail_id'] = main_id
        df_none = df.fillna("")
        self.create_db_table('circrna_identify_detail', df_none.to_dict('r'))
        self.update_db_record('circrna_identify', main_id, insert_dict={'main_id': main_id, 'status': 'end'})
        # if 'circbase' in df_columns:
        #     for i in range(len(df)):
        #         if ',' in df['host_gene_id'][i]:
        #             host_id = df['host_gene_id'][i].split(',')
        #             for l in host_id:
        #                 series = pd.Series({'circrna_id': df['circrna_id'][i], 'host_gene_id': l, 'chr': df['chr'][i],
        #                                     'strand': df['strand'][i],
        #                                     'circrna_start': df['circrna_start'][i],
        #                                     'circrna_end': df['circrna_end'][i],
        #                                     'signal': df['signal'][i], 'circrna_type': df['circrna_type'][i],
        #                                     'circBase': df['circBase'][i]})
        #                 series.name = 'tail_{}'.format(i)
        #                 df = df.append(series)
        #             df.drop([i], inplace=True)
        #     df['detail_id'] = main_id
        #     self.create_db_table('circrna_identify_detail', df.to_dict('r'))
        #     self.update_db_record('circrna_identify', main_id, insert_dict={'main_id': main_id, 'status': 'end'})
        # else:
        #     for i in range(len(df)):
        #         if ',' in df['host_gene_id'][i]:
        #             host_id = df['host_gene_id'][i].split(',')
        #             for l in host_id:
        #                 series = pd.Series({'circrna_id': df['circrna_id'][i], 'host_gene_id': l, 'chr': df['chr'][i],
        #                                     'strand': df['strand'][i],
        #                                     'circrna_start': df['circrna_start'][i],
        #                                     'circrna_end': df['circrna_end'][i],
        #                                     'signal': df['signal'][i], 'circrna_type': df['circrna_type'][i]
        #                                     })
        #                 series.name = 'tail_{}'.format(i)
        #                 df = df.append(series)
        #             df.drop([i], inplace=True)
        #     df['detail_id'] = main_id
        #     self.create_db_table('circrna_identify_detail', df.to_dict('r'))
        #     self.update_db_record('circrna_identify', main_id, insert_dict={'main_id': main_id, 'status': 'end'})


class TestFunction(unittest.TestCase):
    '''
    This is test for the api. Just run this script to do test.
    '''

    def test(self):
        import random
        from mbio.workflows.whole_transcriptome.whole_transcriptome_test_api import WholeTranscriptomeTestApiWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'circrna_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'workflow',
            'name': 'whole_transcriptome.whole_transcriptome_test_api',
            'options': {}
        }
        wheet = Sheet(data=data)
        wf = WholeTranscriptomeTestApiWorkflow(wheet)
        wf.IMPORT_REPORT_DATA = True
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.test_api = wf.api.api('whole_transcriptome.circrna')

        wf.test_api.add_circrna(
            task_id='whole_transcriptome',
            project_sn='whole_transcriptome',
            circ_detail='/mnt/ilustre/users/sanger-dev/workspace/20191219/Circrna_tsg_36603/output/circ_brush/detail.txt',
            taxonomy='Animal'
        )


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)

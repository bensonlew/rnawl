# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import datetime
import json
import unittest

import pandas as pd
from biocluster.api.database.base import report_check

from mbio.api.database.whole_transcriptome.api_base import ApiBase


class LncrnaFamily(ApiBase):
    '''
    last_modify: 2019.11.11
    '''

    def __init__(self, bind_object):
        super(LncrnaFamily, self).__init__(bind_object)

    @report_check
    def add_lncrna_family(self, tabular, params=None):
        task_id = self.bind_object.sheet.id
        project_sn = self.bind_object.sheet.project_sn
        time_now = datetime.datetime.now()
        created_ts = time_now.strftime('%Y-%m-%d %H:%M:%S')
        name = 'Lncrna_family_{}'.format(time_now.strftime('%Y%m%d_%H%M%S'))
        desc = 'Workflow result'
        if not params:
            params = json.dumps({
                'task_id': task_id,
                'submit_location': 'lncrna_family',
                'task_type': 2,
            })
        main_info = {
            'task_id': task_id,
            'project_sn': project_sn,
            'created_ts': created_ts,
            'name': name,
            'desc': desc,
            'params': params,
            'status': 'start',
            'version': 'v1'
        }
        main_id = self.create_db_table('lncrna_family', [main_info])
        df = pd.read_table(tabular)
        df['lnc_fam_id'] = main_id
        self.create_db_table('lncrna_family_detail', df.to_dict('r'))
        self.update_db_record('lncrna_family', main_id, main_id=main_id, status='end')


class TestFunction(unittest.TestCase):
    '''
    This is test for the api. Just run this script to do test.
    '''

    def test(self):
        import random
        from mbio.workflows.lnc_rna.lnc_rna_test_api import LncRnaTestApiWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'lncrna_family_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
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
        wf.test_api = wf.api.api('lnc_rna.lncrna_family')
        wf.test_api.add_lncrna_family(
            '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/lnc_rna/lncrna_family/output/lncRNA_family.tabular'
        )


if __name__ == '__main__':
    unittest.main()

# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

from mbio.api.database.lnc_rna.api_base import ApiBase
from biocluster.api.database.base import report_check
import pandas as pd
from bson.objectid import ObjectId
import unittest

class LncrnaSeqCons(ApiBase):
    '''
    last_modify: 2019.04.02
    '''
    def __init__(self, bind_object):
        super(LncrnaSeqCons, self).__init__(bind_object)

    @report_check
    def add_lncrna_seq_cons(self, tabular, main_id):
        df = pd.read_table(tabular)
        main_id = ObjectId(main_id)
        df['seq_cons_id'] = main_id
        self.create_db_table('sg_lncrna_seq_cons_detail', df.to_dict('r'))
        self.update_db_record('sg_lncrna_seq_cons', main_id, main_id=main_id, status='end')

class TestFunction(unittest.TestCase):
    '''
    This is test for the api. Just run this script to do test.
    '''
    def test(self):
        import random
        from mbio.workflows.lnc_rna.lnc_rna_test_api import LncRnaTestApiWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'lncrna_seq_cons_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
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
        wf.test_api = wf.api.api('lnc_rna.lncrna_seq_cons')
        wf.test_api.add_lncrna_seq_cons(
            tabular='/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/lnc_rna/lncrna_seq_cons/output/lncRNA_sequence_conservation.tabular',
            main_id='5ca2bdeb17b2bf16306027c1'
        )

if __name__ == '__main__':
    unittest.main()
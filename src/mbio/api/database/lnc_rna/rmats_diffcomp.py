# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

from mbio.api.database.lnc_rna.api_base import ApiBase
from biocluster.api.database.base import report_check
import pandas as pd
from bson.objectid import ObjectId
import unittest

class RmatsDiffcomp(ApiBase):
    '''
    last_modify: 2019.03.20
    '''
    def __init__(self, bind_object):
        super(RmatsDiffcomp, self).__init__(bind_object)

    @report_check
    def add_rmats_diffcomp_detail(self, diffcomp_txt, main_id):
        df = pd.read_table(diffcomp_txt)
        df.rna_type.fillna('other', inplace=True)
        df = df.fillna('')
        df['diffcomp_id'] = ObjectId(main_id)
        dict_list = df.to_dict('r')
        for n, d in enumerate(dict_list):
            for k, v in d.items():
                if k.endswith('_JC') or k.endswith('_JCEC'):
                    try:
                        d[k] = float(v)
                    except:
                        pass
            else:
                dict_list[n] = d
        self.create_db_table('sg_splicing_rmats_diffcomp_detail', dict_list)
        insert_dict = {
            'main_id': ObjectId(main_id),
            'status': 'end'
        }
        self.update_db_record('sg_splicing_rmats_diffcomp', main_id, insert_dict=insert_dict)

class TestFunction(unittest.TestCase):
    '''
    This is test for the api. Just run this script to do test.
    '''
    def test(self):
        import random
        from mbio.workflows.lnc_rna.lnc_rna_test_api import LncRnaTestApiWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'rmats_model_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
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
        wf.test_api.add_rmats_model(
            output_dir='/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/lnc_rna/rmats_model/output',
            s3_output=self._sheet.output,
            main_id=''
        )

if __name__ == '__main__':
    unittest.main()
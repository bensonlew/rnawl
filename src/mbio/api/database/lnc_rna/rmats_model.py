# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

from mbio.api.database.lnc_rna.api_base import ApiBase
from biocluster.api.database.base import report_check
import os
from bson.objectid import ObjectId
import unittest

class RmatsModel(ApiBase):
    def __init__(self, bind_object):
        super(RmatsModel, self).__init__(bind_object)

    @report_check
    def add_rmats_model(self, output_dir, s3_output, main_id):
        insert_dict = {
            'pdf_files': [i for i in os.listdir(output_dir) if i.endswith('.pdf')],
            'png_files': [i for i in os.listdir(output_dir) if i.endswith('.png')],
            'graph_dir': s3_output,
            'main_id': ObjectId(main_id),
            'status': 'end'
        }
        self.update_db_record('sg_splicing_rmats_model', main_id, insert_dict=insert_dict)

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
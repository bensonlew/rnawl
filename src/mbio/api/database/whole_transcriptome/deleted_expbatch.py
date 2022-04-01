# -*- coding: utf-8 -*-
# __author__ = 'zoujiaxun'
from bson.objectid import ObjectId
import datetime
import json
import pickle
import unittest
import pandas as pd
from mbio.api.database.whole_transcriptome.api_base import ApiBase
from biocluster.config import Config
import re
import types

class DeletedExpbatch(ApiBase):
    def __init__(self, bind_object):
        super(DeletedExpbatch, self).__init__(bind_object)
        self._project_type = 'whole_transcriptome'

    def deleted_batch(self, task_id, library, level, main_id=None, project_sn='deleted'):
        if main_id is None:
            name = "Deleted_ExpBatch" + '_' + library + '_' + level + '_'
            time_now = datetime.datetime.now()
            name += time_now.strftime("%Y%m%d_%H%M%S")
            main_info = dict(
                project_sn=project_sn,
                task_id=task_id,
                version="v1.1",
                name=name,
                created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
                desc='Deleted ExpBatch Table',
                params=[],
                status="start",
            )
            main_id = self.create_db_table('sg_deleted_expbatch', [main_info])
        else:
            if type(main_id) == str or type(main_id) == bytes or type(main_id) == unicode:
                main_id = ObjectId(main_id)
        record_task = self.db['task'].find_one({'task_id': task_id})
        rnas = [rna for rna, lib in record_task['rna'].items() if lib == library]
        try:
            record_exp_batch = self.db['exp'].find_one({'task_id': task_id, 'level': level, 'is_rmbe': 'true', 'status': 'end'})
            main_exp_id = record_exp_batch['main_id']
            main_exp_id = ObjectId(main_exp_id)
            old_library = list(set(record_exp_batch['library'].split(';')))
            old_category = list(set(record_exp_batch['category'].split(';')))
            new_library = [l for l in old_library if l != library]
            new_library_str = ';'.join(new_library)
            new_category = [c for c in old_category if c not in rnas]
            if not new_category:
                self.update_db_record('exp', main_exp_id, status='failed')
            else:
                new_category_str = ';'.join(new_category)
                self.update_db_record('exp', main_exp_id, library=new_library_str, category=new_category_str)
                self.update_db_record('sg_deleted_expbatch', main_id, library=new_library_str, category=new_category_str)

        except:
            self.bind_object.set_error("删除批次表失败")

class TestFunction(unittest.TestCase):
    '''
    This is test for the api. Just run this script to do test.
    '''
    def test(self):
        import random
        from mbio.workflows.ref_rna_v2.refrna_test_api import RefrnaTestApiWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'add_annotation_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
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
        wf.test_api = wf.api.api('ref_rna_v3.add_annotation')
        wf.test_api.ghost_add_annotation(
            # annotation_file='/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/ref_rna_v3/add_annotation/annotation',
            task_id="tsg_37259",
            # change_database="custom",
            # annot_type='replace',
            # annot_method='all'
        )

if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
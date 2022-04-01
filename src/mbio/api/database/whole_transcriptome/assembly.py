# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import datetime
import json
import pickle
import unittest

from mbio.api.database.whole_transcriptome.api_base import ApiBase


class Assembly(ApiBase):
    def __init__(self, bind_object):
        super(Assembly, self).__init__(bind_object)

    def add_assembly(self, map_dict, task_id, project_sn):
        time_now = datetime.datetime.now()
        name = 'Assembly_{}'.format(time_now.strftime('%Y%m%d_%H%M%S'))
        params = json.dumps({'task_id': task_id, 'submit_location': 'assembly', 'task_type': 2}, sort_keys=True)
        main_dict = {
            'task_id': task_id,
            'project_sn': project_sn,
            'name': name,
            'created_ts': time_now.strftime('%Y-%m-%d %H:%M:%S'),
            'desc': 'Workflow result',
            'params': params,
            'status': 'start',
            'version': 'v1'
        }
        main_id = self.create_db_table('assembly', [main_dict])
        step = self.add_assembly_step(map_dict['step'], main_id)
        code = self.add_assembly_code(map_dict['code'], main_id)
        self.update_db_record('assembly', main_id, insert_dict={
            'main_id': main_id, 'status': 'end', 'step': step, 'code': code
        })

    def add_assembly_step(self, step_file, assembly_id):
        self.bind_object.logger.info('start inserting documents into mongo from {}'.format(step_file))
        docs = pickle.load(open(step_file))
        self.create_db_table('assembly_step', docs, {'assembly_id': assembly_id})
        return list(set(doc['step'] for doc in docs))

    def add_assembly_code(self, code_file, assembly_id):
        self.bind_object.logger.info('start inserting documents into mongo from {}'.format(code_file))
        docs = pickle.load(open(code_file))
        self.create_db_table('assembly_code', docs, {'assembly_id': assembly_id})
        return list(set(doc['class_code'] for doc in docs))


class TestFunction(unittest.TestCase):
    '''
    This is test for the api. Just run this script to do test.
    '''

    def test(self):
        import random
        from mbio.workflows.whole_transcriptome.whole_transcriptome_test_api import WholeTranscriptomeTestApiWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'assembly_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'workflow',
            'name': 'whole_transcriptome.whole_transcriptome_test_api',
            'options': {}
        }
        wheet = Sheet(data=data)
        wf = WholeTranscriptomeTestApiWorkflow(wheet)
        wf.IMPORT_REPORT_DATA = True
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.test_api = wf.api.api('whole_transcriptome.assembly')
        map_dict = {
            'step': '/mnt/ilustre/users/sanger-dev/workspace/20191014/Longrna_workflow_9038_1820/Assembly/output/step.pk',
            'code': '/mnt/ilustre/users/sanger-dev/workspace/20191014/Longrna_workflow_9038_1820/Assembly/output/code.pk'
        }
        wf.test_api.add_assembly(
            map_dict=map_dict,
            task_id='whole_transcriptome',
            project_sn='whole_transcriptome'
        )


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)

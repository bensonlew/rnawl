# -*- coding: utf-8 -*-
# __author__ = 'shicaiping,qinjincheng'

from mainapp.controllers.project.medical_transcriptome_controller import MedicalTranscriptomeController
from mainapp.libs.signature import check_sig
import web
import datetime
import json
import unittest
import os

class RmatsModelAction(MedicalTranscriptomeController):
    '''
    last_modify: 2019.06.18
    '''
    def __init__(self):
        super(RmatsModelAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        args = ['task_id', 'submit_location', 'task_type', 'splicing_id', 'event_type', 'as_type', 'gene_id']
        for arg in args:
            if not hasattr(data, arg):
                info = {'success': False, 'info': 'lack argument -> %s', "variables":[arg], "code" : "C3400402"}
                return json.dumps(info)

        self.medical_transcriptome.delete_main_table('sg_splicing_rmats_model', data.task_id)
        task_id = data.task_id
        task_info = self.medical_transcriptome.get_task_info(task_id=task_id)
        project_sn = task_info['project_sn']
        time_now = datetime.datetime.now()
        rmats_info = self.medical_transcriptome.get_main_info(data.splicing_id, 'sg_splicing_rmats', task_id)
        compare_plan = rmats_info['compare_plan']
        ctrl, test = compare_plan.split('|')
        name = 'Splicing_model_{}_vs_{}_{}'.format(ctrl, test, time_now.strftime('%Y%m%d_%H%M%S'))
        params = json.dumps({
            'task_id': task_id,
            'submit_location': data.submit_location,
            'task_type': int(data.task_type),
            'splicing_id': str(data.splicing_id),
            'event_type': str(data.event_type),
            'as_type': str(data.as_type),
            'gene_id': str(data.gene_id),
        }, sort_keys=True, separators=(',', ':'))
        main_info = {
            'task_id': task_id,
            'name': name,
            'version': "v1",
            'desc': 'alternative splicing model main table',
            'created_ts': time_now.strftime('%Y-%m-%d %H:%M:%S'),
            'params': params,
            'status': 'start',
        }
        main_id = self.medical_transcriptome.insert_main_table('sg_splicing_rmats_model', main_info)
        new_task_id = self.medical_transcriptome.get_new_id(data.task_id)
        main_table_data = {'run_id': new_task_id}

        group_dict = json.dumps(json.loads(rmats_info['params'])['group_dict'])
        result_dir = '{}/'.format(rmats_info['rmats_output'])
        options = {
            's3_file_list': task_id,
            'group_table': group_dict,
            'control_table': compare_plan,
            'result_dir': result_dir,
            'gene_id': data.gene_id,
            'event_type': data.event_type,
            'as_type': data.as_type,
            'main_id': str(main_id),
            'main_table_data': main_table_data,
            'update_info': json.dumps({str(main_id): 'sg_splicing_rmats_model'})
        }
        to_files = [
            'medical_transcriptome_new.medical_transcriptome.export_bam_list(s3_file_list)',
            'medical_transcriptome_new.medical_transcriptome.export_rmats_group_table(group_table)',
            'medical_transcriptome_new.medical_transcriptome.export_rmats_control_table(control_table)'
        ]
        self.set_sheet_data(
            name='medical_transcriptome.report.rmats_model',
            options=options,
            main_table_name=name,
            module_type='workflow',
            to_file=to_files,
            project_sn=project_sn,
            new_task_id=new_task_id,
            task_id=task_id
        )

        run_info = super(RmatsModelAction, self).POST()
        run_info['content'] = {'ids': {'id': str(main_id), 'name': name}}
        return json.dumps(run_info)

class TestFunction(unittest.TestCase):
    '''
    This is test for the controller. Just run this script to do test.
    '''
    def test12(self):
        cmd = 'python /mnt/ilustre/users/sanger-dev/biocluster/bin/webapitest.py post'
        cmd += ' -fr no'
        cmd += ' -c client03'
        cmd += ' -b http://192.168.12.101:9090'
        cmd += ' s/ref_rna_v3/rmats_model'
        args = {
            'task_id': 'tsg_33555',
            'submit_location': 'splicingrmats_model',
            'task_type': '2',
            'splicing_id': '5d033cc217b2bf4c6b88eeda',
            'event_type': 'SE',
            'as_type': 'JC',
            'gene_id': 'ENSG00000166579'
        }
        cmd += ' -n "{}" -d "{}"'.format(';'.join(args.keys()), ';'.join(args.values()))
        print cmd
        os.system(cmd)

if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test12')])
    unittest.TextTestRunner(verbosity=2).run(suite)

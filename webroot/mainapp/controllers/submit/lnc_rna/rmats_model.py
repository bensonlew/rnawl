# -*- coding: utf-8 -*-
# __author__ = 'shicaiping, qinjincheng'

from mainapp.controllers.project.lnc_rna_controller import LncRnaController
from mainapp.libs.signature import check_sig
import web
import datetime
import json
import unittest
import os

class RmatsModelAction(LncRnaController):
    '''
    last_modify: 2019.03.18
    '''
    def __init__(self):
        super(RmatsModelAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        args = ['task_id', 'submit_location', 'task_type', 'splicing_id', 'event_type', 'as_type', 'gene_id']
        for arg in args:
            if not hasattr(data, arg):
                info = {'success': False, 'info': 'lack argument -> {}'.format(arg)}
                return json.dumps(info)

        self.lnc_rna.delete_main_table('sg_splicing_rmats_model', data.task_id)
        task_id = data.task_id
        task_info = self.lnc_rna.get_task_info(task_id=task_id)
        project_sn = task_info['project_sn']
        time_now = datetime.datetime.now()
        created_ts = time_now.strftime('%Y-%m-%d %H:%M:%S')
        rmats_info = self.lnc_rna.get_main_info(data.splicing_id, 'sg_splicing_rmats', task_id)
        compare_plan = rmats_info['compare_plan']
        test, ctrl = compare_plan.split('|')
        name = 'Splicing_model_{}_vs_{}_{}'.format(test, ctrl, time_now.strftime('%Y%m%d_%H%M%S'))
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
            'desc': 'alternative splicing model main table',
            'created_ts': time_now.strftime('%Y-%m-%d %H:%M:%S'),
            'params': params,
            'status': 'start',
        }
        main_id = self.lnc_rna.insert_main_table('sg_splicing_rmats_model', main_info)

        group_dict = json.loads(rmats_info['params'])['group_dict']
        result_dir = '{}/'.format(rmats_info['result_dir'])
        options = {
            's3_file_list': task_id,
            'group_table': group_dict,
            'control_table': compare_plan,
            'result_dir': result_dir,
            'gene_id': data.gene_id,
            'event_type': data.event_type,
            'as_type': data.as_type,
            'main_id': str(main_id),
            'update_info': json.dumps({str(main_id): 'sg_splicing_rmats_model'})
        }
        to_files = [
            'lnc_rna.export_rmats_group_table(group_table)',
            'lnc_rna.export_rmats_control_table(control_table)',
            'lnc_rna.export_bam_list(s3_file_list)'
        ]
        self.set_sheet_data(
            name='lnc_rna.report.rmats_model',
            options=options,
            main_table_name=name,
            module_type='workflow',
            to_file=to_files,
            project_sn=project_sn,
            task_id=task_id
        )

        run_info = super(RmatsModelAction, self).POST()
        run_info['content'] = {'ids': {'id': str(main_id), 'name': name}}
        return json.dumps(run_info)

class TestFunction(unittest.TestCase):
    '''
    This is test for the controller. Just run this script to do test.
    '''
    def test(self):
        cmd = 'python /mnt/ilustre/users/sanger-dev/biocluster/bin/webapitest.py post'
        cmd += ' -fr no'
        cmd += ' -c client03'
        cmd += ' -b http://192.168.12.101:9090'
        cmd += ' s/lnc_rna/rmats_model'
        args = {
            'task_id': 'lnc_rna',
            'submit_location': 'splicingrmats_model',
            'task_type': '2',
            'splicing_id': '5c904b7417b2bf2b0bb46234',
            'event_type': 'SE',
            'as_type': 'JC',
            'gene_id': 'ENSG00000111912',
        }
        cmd += ' -n "{}" -d "{}"'.format(';'.join(args.keys()), ';'.join(args.values()))
        print cmd
        os.system(cmd)

if __name__ == '__main__':
    unittest.main()

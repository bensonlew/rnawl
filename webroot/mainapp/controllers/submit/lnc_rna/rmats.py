# -*- coding: utf-8 -*-
# __author__ = 'shicaiping, qinjincheng'

from mainapp.controllers.project.lnc_rna_controller import LncRnaController
from mainapp.libs.signature import check_sig
import web
import datetime
import json
import unittest
import os

class RmatsAction(LncRnaController):
    '''
    last_modify: 2019.03.19
    '''
    def __init__(self):
        super(RmatsAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        args = ['task_id', 'submit_location', 'task_type', 'group_id', 'group_dict', 'compare_plan', 'control_id']
        for arg in args:
            if not hasattr(data, arg):
                info = {'success': False, 'info': 'lack argument -> {}'.format(arg)}
                return json.dumps(info)

        task_id = data.task_id
        task_info = self.lnc_rna.get_task_info(task_id=task_id)
        project_sn = task_info['project_sn']
        time_now = datetime.datetime.now()
        created_ts = time_now.strftime('%Y-%m-%d %H:%M:%S')
        test, ctrl = data.compare_plan.split('|')
        name = 'Splicing_{}_vs_{}_{}'.format(test, ctrl, time_now.strftime('%Y%m%d_%H%M%S'))
        params = json.dumps({
            'task_id': task_id,
            'submit_location': data.submit_location,
            'task_type': int(data.task_type),
            'group_id': str(data.group_id),
            'group_dict': json.loads(data.group_dict),
            'compare_plan': str(data.compare_plan),
            'control_id': str(data.control_id),
        }, sort_keys=True, separators=(',', ':'))
        main_info = {
            'task_id': task_id,
            'project_sn': project_sn,
            'name': name,
            'desc': 'alternative splicing main table',
            'created_ts': time_now.strftime('%Y-%m-%d %H:%M:%S'),
            'params': params,
            'status': 'start',
            'compare_plan': data.compare_plan
        }
        main_id = self.lnc_rna.insert_main_table('sg_splicing_rmats', main_info)

        ref_gtf = task_info['ref_gtf']
        read_type = task_info['read_type']
        lib_type = 'fr-{}'.format(task_info['lib_type'])
        type_file = task_info['all_gene_type']
        options = {
            'group_table': data.group_dict,
            'control_table': data.compare_plan,
            's3_file_list': task_id,
            'ref_gtf': ref_gtf,
            'read_type': read_type,
            'lib_type': lib_type,
            'splicing_id': str(main_id),
            'type_file': type_file,
            'update_info': json.dumps({str(main_id): 'sg_splicing_rmats'})
        }
        to_files = [
            'lnc_rna.export_rmats_group_table(group_table)',
            'lnc_rna.export_rmats_control_table(control_table)',
            'lnc_rna.export_bam_list(s3_file_list)'
        ]
        self.set_sheet_data(
            name='lnc_rna.report.rmats',
            options=options,
            main_table_name=name,
            module_type='workflow',
            to_file=to_files,
            project_sn=project_sn,
            task_id=task_id
        )

        run_info = super(RmatsAction, self).POST()
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
        cmd += ' s/lnc_rna/rmats'
        args = {
            'task_id': 'tsg_33915',
            'submit_location': 'splicingrmats',
            'task_type': '2',
            'group_id': '5cbd4ccb17b2bf5efca74408',
            'group_dict': json.dumps({'Con':['Con1','Con2', 'Con3', 'Con4', 'Con5'], 'Vit':['Vit1','Vit2', 'Vit3', 'Vit4', 'Vit5']}).replace('"', '\\"'),
            'compare_plan': 'Con|Vit',
            'control_id': '5cbd4ccb17b2bf5efca74409',
        }
        cmd += ' -n "{}" -d "{}"'.format(';'.join(args.keys()), ';'.join(args.values()))
        print cmd
        os.system(cmd)

    def test_3(self):
        cmd = 'python /mnt/ilustre/users/sanger-dev/biocluster/bin/webapitest.py post'
        cmd += ' -fr no'
        cmd += ' -c client03'
        cmd += ' -b http://192.168.12.101:9090'
        cmd += ' s/lnc_rna/rmats'
        args = {
            'task_id': 'lnc_rna',
            'submit_location': 'splicingrmats',
            'task_type': '2',
            'group_id': '5c8f5f4517b2bf10380f1b1a',
            'group_dict': json.dumps({'Con_3':['Con1', 'Con3', 'Con5'], 'Vit_3':['Vit1', 'Vit3', 'Vit5']}).replace('"', '\\"'),
            'compare_plan': 'Con_3|Vit_3',
            'control_id': '5c8f5f4517b2bf10380f1b1b',
        }
        cmd += ' -n "{}" -d "{}"'.format(';'.join(args.keys()), ';'.join(args.values()))
        print cmd
        os.system(cmd)

if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)

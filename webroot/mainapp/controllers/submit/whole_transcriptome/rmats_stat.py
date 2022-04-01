# -*- coding: utf-8 -*-
# __author__ = 'shicaiping, qinjincheng'

import datetime
import json
import os
import unittest

import web

from mainapp.controllers.project.whole_transcriptome_controller import WholeTranscriptomeController
from mainapp.libs.signature import check_sig


class RmatsStatAction(WholeTranscriptomeController):
    '''
    last_modify: 2019.06.14
    '''
    def __init__(self):
        super(RmatsStatAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        args = ['task_id', 'submit_location', 'task_type', 'splicing_id', 'psi', 'pvalue_fdr', 'fdr']
        for arg in args:
            if not hasattr(data, arg):
                info = {'success': False, 'info': 'lack argument -> %s', "variables":[arg], "code" : "C3400502"}
                return json.dumps(info)

        task_id = data.task_id
        task_info = self.whole_transcriptome.get_task_info(task_id=task_id)
        project_sn = task_info['project_sn']
        time_now = datetime.datetime.now()
        ctrl, test = self.whole_transcriptome.get_main_info(
            data.splicing_id, 'splicing_rmats', task_id
        )['compare_plan'].split('|')
        group = {ctrl: 's1', test: 's2'}
        name = 'Splicing_stat_{}_vs_{}_{}'.format(ctrl, test, time_now.strftime('%Y%m%d_%H%M%S'))
        params = json.dumps({
            'task_id': task_id,
            'submit_location': data.submit_location,
            'task_type': int(data.task_type),
            'splicing_id': str(data.splicing_id),
            'psi': data.psi,
            'pvalue_fdr': data.pvalue_fdr,
            'fdr': data.fdr
        }, sort_keys=True, separators=(',', ':'))
        main_info = {
            'task_id': task_id,
            'project_sn': project_sn,
            'name': name,
            'desc': 'alternative splicing stats main table',
            'created_ts': time_now.strftime('%Y-%m-%d %H:%M:%S'),
            'params': params,
            'status': 'start',
            'group': group,
            'version': 'v1'
        }
        main_id = self.whole_transcriptome.insert_main_table('splicing_rmats_stats', main_info)

        options = {
            'root': data.splicing_id,
            'pvalue_fdr': data.pvalue_fdr,
            'fdr': data.fdr,
            'psi': data.psi,
            'main_id': str(main_id),
            'update_info': json.dumps({str(main_id): 'splicing_rmats_stats'})
        }
        to_files = [
            'whole_transcriptome.rmats.export_rmats_root(root)'
        ]
        self.set_sheet_data(
            name='whole_transcriptome.report.rmats_stat',
            options=options,
            main_table_name=name,
            module_type='workflow',
            to_file=to_files,
            project_sn=project_sn,
            task_id=task_id
        )

        run_info = super(RmatsStatAction, self).POST()
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
        cmd += ' s/ref_rna_v3/rmats_stat'
        args = {
            'task_id': 't33555',
            'submit_location': 'splicingrmats_stat',
            'task_type': '2',
            'splicing_id': '5d033cc217b2bf4c6b88eeda',
            'psi': '0.1',
            'pvalue_fdr': 'fdr',
            'fdr': '0.05'
        }
        cmd += ' -n "{}" -d "{}"'.format(';'.join(args.keys()), ';'.join(args.values()))
        print cmd
        os.system(cmd)

    def test34(self):
        cmd = 'python /mnt/ilustre/users/sanger-dev/biocluster/bin/webapitest.py post'
        cmd += ' -fr no'
        cmd += ' -c client03'
        cmd += ' -b http://192.168.12.101:9090'
        cmd += ' s/ref_rna_v3/rmats_stat'
        args = {
            'task_id': 't33555',
            'submit_location': 'splicingrmats_stat',
            'task_type': '2',
            'splicing_id': '5d033cc617b2bf4c6b88eedb',
            'psi': '0.0',
            'pvalue_fdr': 'pvalue',
            'fdr': '0.01'
        }
        cmd += ' -n "{}" -d "{}"'.format(';'.join(args.keys()), ';'.join(args.values()))
        print cmd
        os.system(cmd)

if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test12'), TestFunction('test34')])
    unittest.TextTestRunner(verbosity=2).run(suite)

# -*- coding: utf-8 -*-
# __author__ = 'shicaiping, qinjincheng'

import datetime
import json
import os
import unittest

import web

from mainapp.controllers.project.whole_transcriptome_controller import WholeTranscriptomeController
from mainapp.libs.signature import check_sig


class RmatsDiffcompAction(WholeTranscriptomeController):
    '''
    last_modify: 2019.07.11
    '''
    def __init__(self):
        super(RmatsDiffcompAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        args = ['task_id', 'submit_location', 'task_type', 'main_id',
                'delta_PSI', 'significant_diff', 'significant_value']
        for arg in args:
            if not hasattr(data, arg):
                info = {'success': False, 'info': 'lack argument -> %s', "variables":[arg], "code" : "C3400302"}
                return json.dumps(info)

        task_id = data.task_id
        task_info = self.whole_transcriptome.get_task_info(task_id=task_id)
        project_sn = task_info['project_sn']
        time_now = datetime.datetime.now()
        created_ts = time_now.strftime('%Y-%m-%d %H:%M:%S')
        compare_plans = list()
        for splicing_id in data.main_id.split(','):
            rmats_info = self.whole_transcriptome.get_main_info(splicing_id, 'splicing_rmats', task_id)
            compare_plan = rmats_info['compare_plan']
            compare_plans.append(compare_plan)
        else:
            diffgroups = ','.join(compare_plans)
        name = 'Splicing_diffcomp_{}'.format(time_now.strftime('%Y%m%d_%H%M%S'))
        params = json.dumps({
            'task_id': task_id,
            'submit_location': data.submit_location,
            'task_type': int(data.task_type),
            'main_id': str(data.main_id),
            'delta_PSI': str(data.delta_PSI),
            'significant_diff': str(data.significant_diff),
            'significant_value': str(data.significant_value),
        }, sort_keys=True, separators=(',', ':'))
        main_info = {
            'task_id': task_id,
            'name': name,
            'desc': 'alternative splicing different comparison main table',
            'created_ts': time_now.strftime('%Y-%m-%d %H:%M:%S'),
            'params': params,
            'compare_plans': compare_plans,
            'diffgroups': diffgroups,
            'status': 'start',
            'version': 'v1'
        }
        main_id = self.whole_transcriptome.insert_main_table('splicing_rmats_diffcomp', main_info)

        options = {
            's3_file_list': data.main_id,
            'delta_psi': float(data.delta_PSI),
            'significant_diff': data.significant_diff.lower(),
            'significant_value': float(data.significant_value),
            'main_id': str(main_id),
            'update_info': json.dumps({str(main_id): 'splicing_rmats_diffcomp'})
        }
        to_files = [
            'whole_transcriptome.rmats.export_rmats_detail_path2base(s3_file_list)'
        ]
        self.set_sheet_data(
            name='whole_transcriptome.report.rmats_diffcomp',
            options=options,
            main_table_name=name,
            module_type='workflow',
            to_file=to_files,
            project_sn=project_sn,
            task_id=task_id
        )

        run_info = super(RmatsDiffcompAction, self).POST()
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
        cmd += ' s/ref_rna_v3/rmats_diffcomp'
        args = {
            'task_id': '5da927cf17b2bf11432147f0',
            'submit_location': 'splicingrmats_diffcomp',
            'task_type': '2',
            'main_id': '5da83a4517b2bf0f4d59502a,5da83a7517b2bf0f4d59fac6',
            'delta_PSI': '0',
            'significant_diff': 'pvalue',
            'significant_value': '0.05'
        }
        cmd += ' -n "{}" -d "{}"'.format(';'.join(args.keys()), ';'.join(args.values()))
        print cmd
        os.system(cmd)

if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)

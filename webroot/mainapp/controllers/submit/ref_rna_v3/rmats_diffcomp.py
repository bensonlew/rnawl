# -*- coding: utf-8 -*-
# __author__ = 'shicaiping, qinjincheng'

from mainapp.controllers.project.ref_rna_v2_controller import RefRnaV2Controller
from mainapp.libs.signature import check_sig
import web
import datetime
import json
import unittest
import os

class RmatsDiffcompAction(RefRnaV2Controller):
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
        task_info = self.ref_rna_v2.get_task_info(task_id=task_id)
        project_sn = task_info['project_sn']
        time_now = datetime.datetime.now()
        created_ts = time_now.strftime('%Y-%m-%d %H:%M:%S')
        compare_plans = list()
        for splicing_id in data.main_id.split(','):
            rmats_info = self.ref_rna_v2.get_main_info(splicing_id, 'sg_splicing_rmats', task_id)
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
        }
        main_id = self.ref_rna_v2.insert_main_table('sg_splicing_rmats_diffcomp', main_info)
        new_task_id = self.ref_rna_v2.get_new_id(task_id)
        main_table_data = {'run_id': new_task_id}

        options = {
            's3_file_list': data.main_id,
            'delta_psi': float(data.delta_PSI),
            'significant_diff': data.significant_diff.lower(),
            'significant_value': float(data.significant_value),
            'main_id': str(main_id),
            'main_table_data': main_table_data,
            'update_info': json.dumps({str(main_id): 'sg_splicing_rmats_diffcomp'})
        }
        to_files = [
            'ref_rna_v3.export_rmats_detail_path2base(s3_file_list)'
        ]
        self.set_sheet_data(
            name='ref_rna_v3.report.rmats_diffcomp',
            options=options,
            main_table_name=name,
            module_type='workflow',
            to_file=to_files,
            project_sn=project_sn,
            new_task_id=new_task_id,
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
            'task_id': 'tsg_33555',
            'submit_location': 'splicingrmats_diffcomp',
            'task_type': '2',
            'main_id': '5d033cc217b2bf4c6b88eeda,5d033cc617b2bf4c6b88eedb',
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

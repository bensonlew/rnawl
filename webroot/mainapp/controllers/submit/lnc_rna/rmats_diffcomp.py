# -*- coding: utf-8 -*-
# __author__ = 'shicaiping, qinjincheng'

from mainapp.controllers.project.lnc_rna_controller import LncRnaController
from mainapp.libs.signature import check_sig
import web
import datetime
import json
import unittest
import os

class RmatsDiffcompAction(LncRnaController):
    '''
    last_modify: 2019.03.19
    '''
    def __init__(self):
        super(RmatsDiffcompAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        args = ['task_id', 'submit_location', 'task_type', 'splicing_ids',
                'delta_psi', 'significant_diff', 'significant_value']
        for arg in args:
            if not hasattr(data, arg):
                info = {'success': False, 'info': 'lack argument -> {}'.format(arg)}
                return json.dumps(info)

        task_id = data.task_id
        task_info = self.lnc_rna.get_task_info(task_id=task_id)
        project_sn = task_info['project_sn']
        time_now = datetime.datetime.now()
        created_ts = time_now.strftime('%Y-%m-%d %H:%M:%S')
        compare_plans = list()
        for splicing_id in data.splicing_ids.split(','):
            rmats_info = self.lnc_rna.get_main_info(splicing_id, 'sg_splicing_rmats', task_id)
            compare_plan = rmats_info['compare_plan']
            compare_plans.append(compare_plan)
        name = 'Splicing_different_comparison_{}'.format(time_now.strftime('%Y%m%d_%H%M%S'))
        params = json.dumps({
            'task_id': task_id,
            'submit_location': data.submit_location,
            'task_type': int(data.task_type),
            'splicing_ids': str(data.splicing_ids),
            'delta_psi': str(data.delta_psi),
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
            'status': 'start',
        }
        main_id = self.lnc_rna.insert_main_table('sg_splicing_rmats_diffcomp', main_info)

        type_file = task_info['all_gene_type']
        options = {
            's3_file_list': data.splicing_ids,
            'delta_psi': float(data.delta_psi),
            'significant_diff': data.significant_diff.lower(),
            'significant_value': float(data.significant_value),
            'type_file': type_file,
            'main_id': str(main_id),
            'update_info': json.dumps({str(main_id): 'sg_splicing_rmats_diffcomp'})
        }
        to_files = [
            'lnc_rna.export_rmats_detail_loc2name(s3_file_list)'
        ]
        self.set_sheet_data(
            name='lnc_rna.report.rmats_diffcomp',
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
        cmd += ' s/lnc_rna/rmats_diffcomp'
        args = {
            'task_id': 'lnc_rna',
            'submit_location': 'splicingrmats_diffcomp',
            'task_type': '2',
            'splicing_ids': '5c904b7417b2bf2b0bb46234,5c904ae317b2bf24df83d139',
            'delta_psi': '0.0',
            'significant_diff': 'fdr',
            'significant_value': '0.05',
        }
        cmd += ' -n "{}" -d "{}"'.format(';'.join(args.keys()), ';'.join(args.values()))
        print cmd
        os.system(cmd)

if __name__ == '__main__':
    unittest.main()

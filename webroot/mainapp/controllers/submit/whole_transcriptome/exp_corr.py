# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import datetime
import json
import os
import unittest

import web

from mainapp.controllers.project.whole_transcriptome_controller import WholeTranscriptomeController
from mainapp.libs.signature import check_sig


class ExpCorrAction(WholeTranscriptomeController):
    def __init__(self):
        super(ExpCorrAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        args = ['task_id', 'submit_location', 'task_type', 'library', 'level', 'kind', 'group_id', 'group_dict',
                'take_mean', 'corr_method', 'take_log', 'clus_method', 'dist_method']
        for arg in args:
            if not hasattr(data, arg):
                info = {'success': False, 'info': 'lack argument -> {}'.format(arg)}
                return json.dumps(info)

        task_id = data.task_id
        library = data.library
        level = data.level
        kind = data.kind
        group_id = data.group_id
        group_dict = json.loads(data.group_dict)
        take_mean = data.take_mean
        corr_method = data.corr_method
        take_log = data.take_log
        clus_method = data.clus_method
        dist_method = data.dist_method

        time_now = datetime.datetime.now()
        name = 'Corr_{}_{}'.format(level, time_now.strftime('%Y%m%d_%H%M%S'))
        project_sn = self.whole_transcriptome.get_task_info(task_id)['project_sn']
        params = json.dumps({
            'task_id': task_id,
            'submit_location': data.submit_location,
            'task_type': int(data.task_type),
            'library': library,
            'level': level,
            'kind': kind,
            'group_id': str(group_id),
            'group_dict': group_dict,
            'take_mean': take_mean,
            'corr_method': corr_method,
            'take_log': take_log,
            'dist_method': dist_method,
            'clus_method': clus_method
        }, sort_keys=True, separators=(',', ':'))
        main_dict = {
            'task_id': task_id,
            'project_sn': project_sn,
            'name': name,
            'created_ts': time_now.strftime('%Y-%m-%d %H:%M:%S'),
            'desc': 'Interaction result',
            'params': params,
            'library': library,
            'status': 'start',
            'version': 'v1'
        }
        main_id = self.whole_transcriptome.create_db_table('exp_corr', [main_dict])

        options = {
            'task_id': task_id,
            'library': library,
            'level': level,
            'kind': kind,
            'group_id': group_id,
            'group_dict': data.group_dict,
            'take_mean': take_mean,
            'corr_method': corr_method,
            'dist_method': dist_method,
            'clus_method': clus_method,
            'take_log': take_log,
            'main_id': str(main_id),
            'update_info': json.dumps({str(main_id): 'exp_corr'})
        }
        self.set_sheet_data(
            name='whole_transcriptome.report.exp_corr',
            options=options,
            main_table_name=name,
            module_type='workflow',
            project_sn=project_sn,
            task_id=task_id
        )

        run_info = super(ExpCorrAction, self).POST()
        run_info['content'] = {'ids': {'id': str(main_id), 'name': name}}

        if 'group_id' in data and str(data.group_id).lower() != 'all':
            _ = self.whole_transcriptome.update_group_is_use(data.task_id, data.group_id)
        if 'control_id' in data:
            _ = self.whole_transcriptome.update_group_compare_is_use(data.task_id, data.control_id)

        return json.dumps(run_info)


class TestFunction(unittest.TestCase):
    '''
    This is test for the controller. Just run this script to do test.
    '''

    def test(self):
        cmd = 'python /mnt/ilustre/users/sanger-dev/biocluster/bin/webapitest.py post'
        cmd += ' -fr no'
        cmd += ' -c client03'
        cmd += ' -b http://bcl.tsg.com'
        cmd += ' s/whole_transcriptome/exp_corr'
        args = {
            'task_id': 'tsg_36088',
            'submit_location': 'expcorr',
            'task_type': '2',
            'library': 'long',
            'level': 'G',
            'kind': 'all',
            'group_id': '5dcb545a17b2bf08b8f25ced',
            'group_dict': json.dumps({'Acute': ['A1', 'A2'], 'Stable': ['S1', 'S3'], 'control': ['C2', 'C3']}).replace(
                '"', '\\"'),
            'take_mean': 'no',
            'corr_method': 'kendall',
            'take_log': 'no',
            'clus_method': 'complete',
            'dist_method': 'euclidean',
        }
        cmd += ' -n "{}" -d "{}"'.format(';'.join(args.keys()), ';'.join(args.values()))
        print cmd
        os.system(cmd)


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)

# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import datetime
import json
import os
import unittest

import web

from mainapp.controllers.project.whole_transcriptome_controller import WholeTranscriptomeController
from mainapp.libs.signature import check_sig


class ExpGraphAction(WholeTranscriptomeController):
    def __init__(self):
        super(ExpGraphAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        args = ['task_id', 'submit_location', 'task_type', 'category', 'level', 'kind', 'group_id', 'group_dict']
        for arg in args:
            if not hasattr(data, arg):
                info = {'success': False, 'info': 'lack argument -> {}'.format(arg)}
                return json.dumps(info)

        task_id = data.task_id
        category = data.category
        level = data.level
        kind = data.kind
        group_id = data.group_id
        group_dict = json.loads(data.group_dict)

        time_now = datetime.datetime.now()
        name = 'Graph_{}_{}_{}'.format(level, category, time_now.strftime('%Y%m%d_%H%M%S'))
        project_sn = self.whole_transcriptome.get_task_info(task_id)['project_sn']
        params = json.dumps({
            'task_id': task_id,
            'submit_location': data.submit_location,
            'task_type': int(data.task_type),
            'category': category,
            'level': level,
            'kind': kind,
            'group_id': str(group_id),
            'group_dict': group_dict
        }, sort_keys=True, separators=(',', ':'))
        main_dict = {
            'task_id': task_id,
            'project_sn': project_sn,
            'name': name,
            'created_ts': time_now.strftime('%Y-%m-%d %H:%M:%S'),
            'desc': 'Interaction result',
            'params': params,
            'category': category,
            'status': 'start',
            'version': 'v1'
        }
        main_id = self.whole_transcriptome.create_db_table('exp_graph', [main_dict])

        options = {
            'task_id': task_id,
            'level': level,
            'category': category,
            'kind': kind,
            'group_id': group_id,
            'group_dict': data.group_dict,
            'main_id': str(main_id),
            'update_info': json.dumps({str(main_id): 'exp_graph'})
        }
        self.set_sheet_data(
            name='whole_transcriptome.report.exp_graph',
            options=options,
            main_table_name=name,
            module_type='workflow',
            project_sn=project_sn,
            task_id=task_id
        )

        run_info = super(ExpGraphAction, self).POST()
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
        cmd += ' s/whole_transcriptome/exp_graph'
        args = {
            'task_id': 'tsg_36088',
            'submit_location': 'expgraph',
            'task_type': '2',
            'category': 'mRNA',
            'level': 'T',
            'kind': 'ref',
            'group_id': '5dcb545a17b2bf08b8f25ced',
            'group_dict': json.dumps({'Acute': ['A1', 'A2'], 'Stable': ['S1', 'S3'], 'control': ['C2', 'C3']}).replace(
                '"', '\\"'),
        }
        cmd += ' -n "{}" -d "{}"'.format(';'.join(args.keys()), ';'.join(args.values()))
        print cmd
        os.system(cmd)


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)

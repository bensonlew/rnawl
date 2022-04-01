# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import datetime
import json
import os
import unittest
from bson.objectid import ObjectId
import web

from mainapp.controllers.project.whole_transcriptome_controller import WholeTranscriptomeController
from mainapp.libs.signature import check_sig


class ExpVennAction(WholeTranscriptomeController):
    def __init__(self):
        super(ExpVennAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        args = ['task_id', 'submit_location', 'task_type', 'category', 'level', 'kind', 'group_id', 'group_dict',
                'threshold', 'exp_id']
        for arg in args:
            if not hasattr(data, arg):
                info = {'success': False, 'info': 'lack argument -> {}'.format(arg)}
                return json.dumps(info)

        exp_dict = self.whole_transcriptome.db['exp'].find_one({'task_id': data.task_id, "level": data.level, 'main_id': ObjectId(data.exp_id)})
        is_rmbe = str(exp_dict['is_rmbe']).lower()
        if 'is_rmbe' not in exp_dict or is_rmbe == 'false':
            exp_id = str(exp_dict['main_id'])
        if exp_dict['is_rmbe'] == 'true':
            exp_id = str(exp_dict['batch_main_id'])

        task_id = data.task_id
        category = data.category
        level = data.level
        kind = data.kind
        group_id = data.group_id
        group_dict = json.loads(data.group_dict)
        threshold = float(data.threshold)

        time_now = datetime.datetime.now()
        name = 'Venn_{}_{}_{}'.format(level, category, time_now.strftime('%Y%m%d_%H%M%S'))
        project_sn = self.whole_transcriptome.get_task_info(task_id)['project_sn']
        params = json.dumps({
            'task_id': task_id,
            'submit_location': data.submit_location,
            'task_type': int(data.task_type),
            'category': category,
            'level': level,
            'kind': kind,
            'group_id': str(group_id),
            'group_dict': group_dict,
            'threshold': data.threshold,
            'exp_id': data.exp_id
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
            'version': 'v1.1'
        }
        main_id = self.whole_transcriptome.create_db_table('exp_venn', [main_dict])

        options = {
            'task_id': task_id,
            'level': level,
            'category': category,
            'kind': kind,
            'group_id': group_id,
            'group_dict': data.group_dict,
            'threshold': threshold,
            'main_id': str(main_id),
            'update_info': json.dumps({str(main_id): 'exp_venn'}),
            'is_rmbe': is_rmbe
        }
        self.set_sheet_data(
            name='whole_transcriptome_v1_1.report.exp_venn',
            options=options,
            main_table_name=name,
            module_type='workflow',
            project_sn=project_sn,
            task_id=task_id
        )

        run_info = super(ExpVennAction, self).POST()
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
        cmd += ' s/whole_transcriptome_v1_1/exp_venn'
        args = {
            'task_id': 'ctq8md746f5i66sntrpktbdijm',
            'submit_location': 'expvenn',
            'task_type': '2',
            'category': 'mRNA',
            'level': 'T',
            'kind': 'ref',
            'group_id': '5f7366ec17b2bf6ed1ed37a4',
            'group_dict': json.dumps({'NFD': ['NFD1', 'NFD2', 'NFD3', 'NFD4'], 'HFD': ['HFD1', 'HFD2', 'HFD3', 'HFD4'], 'NAC_HFD': ['NAC_HFD1', 'NAC_HFD2', 'NAC_HFD3', 'NAC_HFD4']}).replace(
                '"', '\\"'),
            'threshold': '3.0',
            'exp_id': '5fa2587317b2bf223dc1ecb7'
        }
        cmd += ' -n "{}" -d "{}"'.format(';'.join(args.keys()), ';'.join(args.values()))
        print cmd
        os.system(cmd)


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)

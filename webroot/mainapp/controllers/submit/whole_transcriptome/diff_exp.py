# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import datetime
import json
import os
import unittest
from collections import OrderedDict

import web

from mainapp.controllers.project.whole_transcriptome_controller import WholeTranscriptomeController
from mainapp.libs.signature import check_sig


class DiffExpAction(WholeTranscriptomeController):
    def __init__(self):
        super(DiffExpAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        args = ['task_id', 'submit_location', 'task_type', 'category', 'level', 'kind', 'filter', 'threshold',
                'group_id', 'group_dict', 'control_id', 'stat_type', 'stat_cutoff', 'fc', 'diff_method']
        for arg in args:
            if not hasattr(data, arg):
                info = {'success': False, 'info': 'lack argument -> {}'.format(arg)}
                return json.dumps(info)

        group_dict = json.loads(data.group_dict, object_pairs_hook=OrderedDict)
        if str(data.group_id).lower() == 'all':
            samples = group_dict['all']
            group_dict = OrderedDict([(x, [x]) for x in samples])
        group_size = list()
        success = list()
        for key in group_dict:
            group_size.append(len(group_dict[key]))
        group_size.sort()
        if group_size[0] == group_size[-1] == 1:
            if data.diff_method == 'DESeq2':
                success.append('比较组别中存在无生物学重复的情况，请选择DEGseq或者edgeR软件运行！')
        elif group_size[0] == 1 and group_size[-1] >= 2:
            if data.diff_method == 'DESeq2':
                success.append('比较组别中存在无生物学重复的情况，请选择DEGseq或者edgeR软件运行！')
        # elif group_size[0] >= 2:
        #     if data.diff_method == 'DEGseq':
        #         success.append('If your experiment has biological replicates, please select DESeq2 or edgeR')
        if success:
            variables = list()
            variables.append(success[0])
            info = {'success': False, 'info': '%s' % success[0], 'code': 'C2900302', 'variables': variables}
            return json.dumps(info)

        task_id = data.task_id
        category = data.category
        level = data.level
        kind = data.kind
        filter = data.filter
        threshold = data.threshold
        group_id = data.group_id
        group_dict = json.loads(data.group_dict)
        control_id = data.control_id
        stat_type = data.stat_type
        stat_cutoff = data.stat_cutoff
        fc = data.fc
        diff_method = data.diff_method

        time_now = datetime.datetime.now()
        name = 'DE_{}_{}_{}'.format(level, category, time_now.strftime('%Y%m%d_%H%M%S'))
        project_sn = self.whole_transcriptome.get_task_info(task_id)['project_sn']
        param_dict = {
            'task_id': task_id,
            'submit_location': data.submit_location,
            'task_type': int(data.task_type),
            'category': category,
            'level': level,
            'kind': kind,
            'filter': filter,
            'threshold': str(threshold),
            'group_id': str(group_id),
            'group_dict': group_dict,
            'control_id': str(control_id),
            'stat_type': stat_type,
            'stat_cutoff': str(stat_cutoff),
            'fc': str(fc),
            'diff_method': diff_method
        }
        if hasattr(data, 'background'):
            param_dict['background'] = data.background
        if hasattr(data, 'correct_method'):
            param_dict['correct_method'] = data.correct_method
        params = json.dumps(param_dict, sort_keys=True, separators=(',', ':'))
        main_dict = {
            'task_id': task_id,
            'project_sn': project_sn,
            'name': name,
            'created_ts': time_now.strftime('%Y-%m-%d %H:%M:%S'),
            'desc': 'Interaction result',
            'params': params,
            'level': level,
            'category': category,
            'status': 'start',
            'version': 'v1'
        }
        main_id = self.whole_transcriptome.create_db_table('diff', [main_dict])

        options = {
            'task_id': task_id,
            'category': category,
            'level': level,
            'kind': kind,
            'group_id': group_id,
            'group_dict': data.group_dict,
            'control_id': control_id,
            'program': diff_method,
            'filter': filter,
            'threshold': threshold,
            'stat_type': stat_type,
            'stat_cutoff': stat_cutoff,
            'fc': fc,
            'main_id': str(main_id),
            'update_info': json.dumps({str(main_id): 'diff'})
        }
        #因页面参数隐藏，取消该参数，所有差异分析均以category指定的rna为背景，不再考虑mRNA或其他mRNA
        if hasattr(data, 'background'):
            options['background'] = data.background
        if hasattr(data, 'correct_method'):
            options['method'] = data.correct_method
        self.set_sheet_data(
            name='whole_transcriptome.report.diff_exp',
            options=options,
            main_table_name=name,
            module_type='workflow',
            project_sn=project_sn,
            task_id=task_id
        )

        run_info = super(DiffExpAction, self).POST()
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
        cmd += ' s/whole_transcriptome/diff_exp'
        args = {
            'task_id': 'tsg_36088',
            'submit_location': 'diffexp',
            'task_type': '2',
            'category': 'mRNA',
            'level': 'G',
            'kind': 'all',
            'background': 'mRNA,lncRNA',
            'filter': 'none',
            'threshold': '0.0',
            'group_id': '5dd1f88b17b2bf10f6f24eb6',
            'group_dict': json.dumps({'Acute': ['A1', 'A2'], 'Stable': ['S1', 'S3'], 'control': ['C2', 'C3']}).replace(
                '"', '\\"'),
            'control_id': '5dd1f88b17b2bf10f6f24eb7',
            'stat_type': 'padjust',
            'stat_cutoff': '0.001',
            'fc': '2.0',
            'diff_method': 'DEGseq',
            'correct_method': 'BH'
        }
        cmd += ' -n "{}" -d "{}"'.format(';'.join(args.keys()), ';'.join(args.values()))
        print cmd
        os.system(cmd)


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)

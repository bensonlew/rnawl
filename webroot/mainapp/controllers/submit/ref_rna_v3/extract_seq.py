# -*- coding: utf-8 -*-
# __author__ = 'shicaiping'

from mainapp.controllers.project.ref_rna_v2_controller import RefRnaV2Controller
from mainapp.libs.signature import check_sig
import web
import datetime
import json
import unittest
import os


class ExtractSeqAction(RefRnaV2Controller):
    '''
    last_modify: 20200707
    '''
    def __init__(self):
        super(ExtractSeqAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        args = ['task_id', 'submit_location', 'task_type', 'extract_info']
        for arg in args:
            if not hasattr(data, arg):
                info = {'success': False, 'info': 'lack argument -> %s', "variables":[arg], "code" : "C3400102"}
                return json.dumps(info)

        task_id = data.task_id
        task_info = self.ref_rna_v2.get_task_info(task_id=task_id)
        project_sn = task_info['project_sn']
        time_now = datetime.datetime.now()
        created_ts = time_now.strftime('%Y-%m-%d %H:%M:%S')
        name = 'ExtractSeq_{}'.format(time_now.strftime('%Y%m%d_%H%M%S'))
        params_dict = dict()
        for each in args:
            if each == "task_type":
                params_dict[each] = int(data[each])
            elif each == 'extract_info':
                params_dict[each] = json.loads(data[each])
            else:
                params_dict[each] = data[each]
        params = json.dumps(params_dict, sort_keys=True, separators=(',', ':'))
        main_info = {
            'task_id': task_id,
            'project_sn': project_sn,
            'name': name,
            'desc': 'extract sequence main table',
            'created_ts': created_ts,
            'params': params,
            'status': 'start',
            'version': 'v3.1'
        }
        main_id = self.ref_rna_v2.insert_main_table('sg_extract', main_info)
        new_task_id = self.ref_rna_v2.get_new_id(task_id)
        main_table_data = {'run_id': new_task_id}
        try:
            task_info['fq_type']
            options = {
                'bam_file_list': task_id,
                'seq_type': task_info['fq_type'],
                'extract_info': data.extract_info,
                'task_id': data.task_id,
                'main_id': str(main_id),
                'main_table_data': main_table_data,
                'update_info': json.dumps({str(main_id): 'sg_extract'})
            }
        except:
            options = {
                'bam_file_list': task_id,
                'extract_info': data.extract_info,
                'task_id': data.task_id,
                'main_id': str(main_id),
                'main_table_data': main_table_data,
                'update_info': json.dumps({str(main_id): 'sg_extract'})
            }
        to_files = [
            'ref_rna_v3.export_bam_list(bam_file_list)'
        ]
        self.set_sheet_data(
            name='ref_rna_v3.report.extract_seq',
            options=options,
            main_table_name=name,
            module_type='workflow',
            to_file=to_files,
            project_sn=project_sn,
            task_id=task_id,
            new_task_id=new_task_id,
        )

        run_info = super(ExtractSeqAction, self).POST()
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
        cmd += ' -b http://bcl.tsg.com'
        cmd += ' s/ref_rna_v3/extract_seq'
        args = {
            'task_id': 'tsg_37858',
            'submit_location': 'extract_seq',
            'task_type': '2',
            'extract_info': json.dumps({"bam": ["W12", "W13", "W14"], "mapped": ["W12", "W13"]}).replace('"', '\\"'),
        }
        cmd += ' -n "{}" -d "{}"'.format(';'.join(args.keys()), ';'.join(args.values()))
        print cmd
        os.system(cmd)


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)


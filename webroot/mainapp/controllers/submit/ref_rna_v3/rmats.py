# -*- coding: utf-8 -*-
# __author__ = 'shicaiping, qinjincheng'

from mainapp.controllers.project.ref_rna_v2_controller import RefRnaV2Controller
from mainapp.libs.signature import check_sig
import web
import datetime
import json
import unittest
import os

class RmatsAction(RefRnaV2Controller):
    '''
    last_modify: 2019.06.14
    '''
    def __init__(self):
        super(RmatsAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        args = ['task_id', 'submit_location', 'task_type', 'group_id', 'group_dict', 'control_id', 'compare_plan']
        for arg in args:
            if not hasattr(data, arg):
                info = {'success': False, 'info': 'lack argument -> %s', "variables":[arg], "code" : "C3400102"}
                return json.dumps(info)

        task_id = data.task_id
        task_info = self.ref_rna_v2.get_task_info(task_id=task_id)
        project_sn = task_info['project_sn']
        time_now = datetime.datetime.now()
        created_ts = time_now.strftime('%Y-%m-%d %H:%M:%S')
        ctrl, test = data.compare_plan.split('|')
        name = 'Splicing_{}_vs_{}_{}'.format(ctrl, test, time_now.strftime('%Y%m%d_%H%M%S'))
        params = json.dumps({
            'task_id': task_id,
            'submit_location': data.submit_location,
            'task_type': int(data.task_type),
            'group_id': str(data.group_id),
            'group_dict': json.loads(data.group_dict),
            'control_id': str(data.control_id),
            'compare_plan': str(data.compare_plan)
        }, sort_keys=True, separators=(',', ':'))
        main_info = {
            'task_id': task_id,
            'project_sn': project_sn,
            'name': name,
            'version':"v3.3",#modify by fwy,20210830,rmats专属版本,为了跟之前版本进行区分
            'desc': 'alternative splicing main table',
            'created_ts': time_now.strftime('%Y-%m-%d %H:%M:%S'),
            'params': params,
            'status': 'start',
            'compare_plan': data.compare_plan,
            "samples": {
                "s1": data.compare_plan.split("|")[-1],
                "s2": data.compare_plan.split("|")[0]
            }
        }
        main_id = self.ref_rna_v2.insert_main_table('sg_splicing_rmats', main_info)
        new_task_id = self.ref_rna_v2.get_new_id(task_id)
        main_table_data = {'run_id': new_task_id}

        as_gtf = task_info['as_gtf']
        if task_info['fq_type'] == 'PE':
            seq_type = 'paired'
        elif task_info['fq_type'] == 'SE':
            seq_type = 'single'
        if task_info['strand_specific']:
            if task_info['strand_dir'][0] == 'R':
                lib_type = 'fr-firststrand'
            elif task_info['strand_dir'][0] == 'F':
                lib_type = 'fr-secondstrand'
        else:
            lib_type = 'fr-unstranded'
        options = {
            's3_file_list': task_id,
            'group_table': data.group_dict,
            'control_table': data.compare_plan,
            'ref_gtf': as_gtf,
            'seq_type': seq_type,
            'lib_type': lib_type,
            'splicing_id': str(main_id),
            'task_id': data.task_id,
            'main_table_data': main_table_data,
            'update_info': json.dumps({str(main_id): 'sg_splicing_rmats'})
        }
        to_files = [
            'ref_rna_v3.export_rmats_group_table(group_table)',
            'ref_rna_v3.export_rmats_control_table(control_table)',
            'ref_rna_v3.export_bam_list(s3_file_list)'
        ]
        self.set_sheet_data(
            name='ref_rna_v3.report.rmats',
            options=options,
            main_table_name=name,
            module_type='workflow',
            to_file=to_files,
            project_sn=project_sn,
            task_id=task_id,
            new_task_id=new_task_id,
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
        cmd += ' s/ref_rna_v3/rmats'
        args = {
            'task_id': 'tsg_33555',
            'submit_location': 'splicingrmats',
            'task_type': '2',
            'group_id': '5d01f6f6c6598df70a8b456d',
            'group_dict': json.dumps({'ctrl':['Con1', 'Con3', 'Con5'], 'test':['Vit1', 'Vit3', 'Vit5']}).replace('"', '\\"'),
            'compare_plan': 'ctrl|test',
            'control_id': '5d01f706c6598dfe0a8b4568',
        }
        cmd += ' -n "{}" -d "{}"'.format(';'.join(args.keys()), ';'.join(args.values()))
        print cmd
        os.system(cmd)

    def test12(self):
        cmd = 'python /mnt/ilustre/users/sanger-dev/biocluster/bin/webapitest.py post'
        cmd += ' -fr no'
        cmd += ' -c client03'
        cmd += ' -b http://192.168.12.101:9090'
        cmd += ' s/ref_rna_v3/rmats'
        args = {
            'task_id': 'tsg_33555',
            'submit_location': 'splicingrmats',
            'task_type': '2',
            'group_id': '5d01b7e3c6598da50b8b456e',
            'group_dict': json.dumps({'Con_12':['Con1', 'Con2'], 'Vit_12':['Vit1', 'Vit2']}).replace('"', '\\"'),
            'compare_plan': 'Con_12|Vit_12',
            'control_id': '5d01da39c6598d540b8b4570',
        }
        cmd += ' -n "{}" -d "{}"'.format(';'.join(args.keys()), ';'.join(args.values()))
        print cmd
        os.system(cmd)

    def test34(self):
        cmd = 'python /mnt/ilustre/users/sanger-dev/biocluster/bin/webapitest.py post'
        cmd += ' -fr no'
        cmd += ' -c client03'
        cmd += ' -b http://192.168.12.101:9090'
        cmd += ' s/ref_rna_v3/rmats'
        args = {
            'task_id': 'tsg_33555',
            'submit_location': 'splicingrmats',
            'task_type': '2',
            'group_id': '5d01b80bc6598d220b8b456d',
            'group_dict': json.dumps({'Con_34':['Con3', 'Con4'], 'Vit_34':['Vit3', 'Vit4']}).replace('"', '\\"'),
            'compare_plan': 'Con_34|Vit_34',
            'control_id': '5d01da52c6598d890b8b456b',
        }
        cmd += ' -n "{}" -d "{}"'.format(';'.join(args.keys()), ';'.join(args.values()))
        print cmd
        os.system(cmd)

if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test12'), TestFunction('test34')])
    unittest.TextTestRunner(verbosity=2).run(suite)

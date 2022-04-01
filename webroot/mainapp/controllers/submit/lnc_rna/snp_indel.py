# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

from mainapp.controllers.project.lnc_rna_controller import LncRnaController
from mainapp.libs.signature import check_sig
import web
import datetime
import json
import unittest
import os

class SnpIndelAction(LncRnaController):
    '''
    last_modify: 2019.04.23
    '''
    def __init__(self):
        super(SnpIndelAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        args = ['task_id', 'submit_location', 'task_type', 'method_type']
        for arg in args:
            if not hasattr(data, arg):
                info = {'success': False, 'info': 'lack argument -> {}'.format(arg)}
                return json.dumps(info)

        task_id = data.task_id
        task_info = self.lnc_rna.get_task_info(task_id=task_id)
        project_sn = task_info['project_sn']
        time_now = datetime.datetime.now()
        created_ts = time_now.strftime('%Y-%m-%d %H:%M:%S')
        name = 'Snp_{}_{}'.format(data.method_type, time_now.strftime('%Y%m%d_%H%M%S'))
        desc = 'Snp main table'
        params = json.dumps({
            'method_type': data.method_type,
            'submit_location': data.submit_location,
            'task_id': task_id,
            'task_type': int(data.task_type),
        }, sort_keys=True, separators=(',', ':'))
        main_info = {
            'task_id': task_id,
            'project_sn': project_sn,
            'created_ts': created_ts,
            'name': name,
            'desc': desc,
            'params': params,
            'status': 'start'
        }
        main_id = self.lnc_rna.insert_main_table('sg_snp', main_info)

        ref_gtf = task_info['ref_gtf']
        ref_fa = task_info['ref_genome']
        ref_genome = task_info['organism_name']
        des = task_info['des']
        des_type = task_info['des_type']
        type_file = task_info['all_gene_type']
        options = {
            'method_type': data.method_type,
            's3_file_list': task_id,
            'ref_fa': ref_fa,
            'ref_gtf': ref_gtf,
            'ref_genome': ref_genome,
            'des': des,
            'des_type': des_type,
            'main_id': str(main_id),
            'type_file': type_file,
            'update_info': json.dumps({str(main_id): 'sg_snp'})
        }
        to_files = ['lnc_rna.export_bam_list(s3_file_list)']
        self.set_sheet_data(
            name='lnc_rna.report.snp_indel',
            options=options,
            main_table_name=name,
            module_type='workflow',
            to_file=to_files,
            project_sn=project_sn,
            task_id=task_id
        )

        run_info = super(SnpIndelAction, self).POST()
        run_info['content'] = {'ids': {'id': str(main_id), 'name': name}}
        return json.dumps(run_info)

class TestFunction(unittest.TestCase):
    '''
    This is test for the controller. Just run this script to do test.
    '''
    def test_samtools(self):
        cmd = 'python /mnt/ilustre/users/sanger-dev/biocluster/bin/webapitest.py post'
        cmd += ' -fr no'
        cmd += ' -c client03'
        cmd += ' -b http://192.168.12.101:9090'
        cmd += ' s/lnc_rna/snp_indel'
        args = {
            'task_id': 'lnc_rna',
            'task_type': '2',
            'submit_location': 'snp',
            'method_type': 'samtools'
        }
        cmd += ' -n "{}" -d "{}"'.format(';'.join(args.keys()), ';'.join(args.values()))
        print cmd
        os.system(cmd)

    def test_gatk(self):
        cmd = 'python /mnt/ilustre/users/sanger-dev/biocluster/bin/webapitest.py post'
        cmd += ' -fr no'
        cmd += ' -c client03'
        cmd += ' -b http://192.168.12.101:9090'
        cmd += ' s/lnc_rna/snp_indel'
        args = {
            'task_id': 'tsg_33915',
            'task_type': '2',
            'submit_location': 'snp',
            'method_type': 'gatk'
        }
        cmd += ' -n "{}" -d "{}"'.format(';'.join(args.keys()), ';'.join(args.values()))
        print cmd
        os.system(cmd)

if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test_gatk')])
    unittest.TextTestRunner(verbosity=2).run(suite)

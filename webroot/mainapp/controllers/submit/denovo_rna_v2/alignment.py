# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import datetime
import json
import os
import unittest

import web

from mainapp.controllers.project.denovo_rna_v2_controller import DenovoRnaV2Controller
from mainapp.libs.signature import check_sig


class AlignmentAction(DenovoRnaV2Controller):
    '''
    last_modify: 2019.10.17
    '''

    def __init__(self):
        super(AlignmentAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        args = ['task_id', 'submit_location', 'task_type', 'query_source', 'subject_source']
        if data.query_source == 'denovo':
            args += ['exp_level', 'geneset_id', 'query_type']
            if data.subject_source == 'ref':
                args += ['subject_target', 'organism_name', 'genome_id']
            elif data.subject_source == 'upload':
                args += ['subject_file', 'subject_id', 'subject_type']
        elif data.query_source == 'upload' and data.subject_source == 'denovo':
            args += ['query_file', 'query_id', 'query_type', 'subject_target']
        args += ['program', 'evalue', 'max_hsps']
        for arg in args:
            if not hasattr(data, arg):
                info = {'success': False, 'info': 'lack argument -> (%s)', "variables":[arg], "code" : "C1602801"}
                return json.dumps(info)

        for document in self.denovo_rna_v2.db['sg_blast'].find({'task_id': data.task_id}):
            _id = document['_id']
            self.denovo_rna_v2.db['sg_blast'].remove({'main_id': _id})
            self.denovo_rna_v2.db['sg_blast_detail'].delete_many({'blast_id': _id})
            self.denovo_rna_v2.db['sg_status'].remove({'task_id': data.task_id, 'table_id': _id})
        else:
            print('succeed in cleaning blast related records by task id ({})'.format(data.task_id))

        task_id = data.task_id
        task_info = self.denovo_rna_v2.get_task_info(task_id)
        project_sn = task_info['project_sn']
        time_now = datetime.datetime.now()
        created_ts = time_now.strftime('%Y-%m-%d %H:%M:%S')
        name = 'Blast_{}_{}_{}'.format(data.query_source, data.subject_source, time_now.strftime('%Y%m%d_%H%M%S'))
        desc = 'Blast main table'
        params = json.dumps(
            dict((arg, int(getattr(data, arg))) if arg == 'task_type' else (arg, getattr(data, arg)) for arg in args),
            sort_keys=True, separators=(',', ':'))
        main_info = {
            'task_id': task_id,
            'project_sn': project_sn,
            'name': name,
            'desc': desc,
            'created_ts': created_ts,
            'params': params,
            'status': 'start'
        }
        main_id = self.denovo_rna_v2.insert_main_table('sg_blast', main_info)

        options = {
            'query_source': data.query_source,
            'subject_source': data.subject_source,
            'program': data.program,
            'evalue': data.evalue,
            'max_hsps': data.max_hsps,
            'main_id': str(main_id),
            'update_info': json.dumps({str(main_id): 'sg_blast'})
        }
        if data.query_source == 'denovo':
            options.update({
                'exp_level': data.exp_level,
                'geneset_id': data.geneset_id,
                'query_type': data.query_type
            })
            if data.subject_source == 'ref':
                options.update({
                    'subject_target': data.subject_target,
                    'subject_type': 'prot' if data.subject_target == 'peptide' else 'nucl',
                    'genome_id': data.genome_id
                })
            elif data.subject_source == 'upload':
                options.update({
                    'subject_file': data.subject_file,
                    'subject_type': data.subject_type
                })
        elif data.query_source == 'upload' and data.subject_source == 'denovo':
            options.update({
                'query_file': data.query_file,
                'query_type': data.query_type,
                'subject_target': data.subject_target,
            })
        self.set_sheet_data(
            name='denovo_rna_v2.report.alignment',
            options=options,
            main_table_name=name,
            module_type='workflow',
            project_sn=project_sn,
            task_id=task_id)

        run_info = super(AlignmentAction, self).POST()
        run_info['content'] = {'ids': {'id': str(main_id), 'name': name}}
        if 'group_id' in data and str(data.group_id).lower() != 'all':
            _ = self.denovo_rna_v2.update_group_is_use(data.task_id, data.group_id)
        if 'control_id' in data:
            _ = self.denovo_rna_v2.update_group_compare_is_use(data.task_id, data.control_id)
        return json.dumps(run_info)


class TestFunction(unittest.TestCase):
    '''
    This is test for the controller. Just run this script to do test.
    '''

    def test_denovo_ref_genome(self):
        cmd = 'python /mnt/ilustre/users/sanger-dev/biocluster/bin/webapitest.py post'
        cmd += ' -fr no'
        cmd += ' -c client03'
        # cmd += ' -b http://192.168.12.101:9090'
        cmd += ' -b http://bcl.tsg.com'
        cmd += ' s/denovo_rna_v2/alignment'
        args = {
            'task_id': 'tsg_33857',
            'submit_location': 'blast',
            'task_type': '2',
            'query_source': 'denovo',
            'exp_level': 'G',
            'geneset_id': '5d37b5f5c6598d2f578b456c',
            'query_type': 'nucl',
            'subject_source': 'ref',
            'subject_target': 'genome',
            'organism_name': 'Saccharomyces_cerevisiae',
            'genome_id': 'GM0265',
            'program': 'blat',
            'evalue': '1e-1',
            'max_hsps': '10',
        }
        cmd += ' -n "{}" -d "{}"'.format(';'.join(args.keys()), ';'.join(args.values()))
        print cmd
        os.system(cmd)

    def test_denovo_ref_database(self):
        cmd = 'python /mnt/ilustre/users/sanger-dev/biocluster/bin/webapitest.py post'
        cmd += ' -fr no'
        cmd += ' -c client03'
        # cmd += ' -b http://192.168.12.101:9090'
        cmd += ' -b http://bcl.tsg.com'
        cmd += ' s/denovo_rna_v2/alignment'
        args = {
            'task_id': 'tsg_33857',
            'submit_location': 'blast',
            'task_type': '2',
            'query_source': 'denovo',
            'exp_level': 'G',
            'geneset_id': '5d37b5f5c6598d2f578b456c',
            'query_type': 'prot',
            'subject_source': 'ref',
            'subject_target': 'peptide',
            'organism_name': 'Saccharomyces_cerevisiae',
            'genome_id': 'GM0265',
            'program': 'blastp',
            'evalue': '1e-1',
            'max_hsps': '10',
        }
        cmd += ' -n "{}" -d "{}"'.format(';'.join(args.keys()), ';'.join(args.values()))
        print cmd
        os.system(cmd)

    def test_denovo_upload(self):
        cmd = 'python /mnt/ilustre/users/sanger-dev/biocluster/bin/webapitest.py post'
        cmd += ' -fr no'
        cmd += ' -c client03'
        # cmd += ' -b http://192.168.12.101:9090'
        cmd += ' -b http://bcl.tsg.com'
        cmd += ' s/denovo_rna_v2/alignment'
        args = {
            'task_id': 'tsg_33857',
            'submit_location': 'blast',
            'task_type': '2',
            'query_source': 'denovo',
            'exp_level': 'G',
            'geneset_id': '5d37b5f5c6598d2f578b456c',
            'query_type': 'prot',
            'subject_source': 'upload',
            # 'upload_method': 'file',
            'subject_file': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/denovo_rna_v2/transcript.fa',
            'subject_id': '123456789',
            'subject_type': 'nucl',
            'program': 'tblastn',
            'evalue': '1e-1',
            'max_hsps': '10',
        }
        cmd += ' -n "{}" -d "{}"'.format(';'.join(args.keys()), ';'.join(args.values()))
        print cmd
        os.system(cmd)

    def test_upload_denovo(self):
        cmd = 'python /mnt/ilustre/users/sanger-dev/biocluster/bin/webapitest.py post'
        cmd += ' -fr no'
        cmd += ' -c client03'
        # cmd += ' -b http://192.168.12.101:9090'
        cmd += ' -b http://bcl.tsg.com'
        cmd += ' s/denovo_rna_v2/alignment'
        args = {
            'task_id': 'tsg_33857',
            'submit_location': 'blast',
            'task_type': '2',
            'query_source': 'upload',
            # 'upload_method': 'file',
            'query_file': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/denovo_rna_v2/cds.fa',
            'query_id': '987654321',
            'query_type': 'nucl',
            'subject_source': 'denovo',
            'subject_target': 'filter',
            'program': 'blastn',
            'evalue': '1e-1',
            'max_hsps': '10',
        }
        cmd += ' -n "{}" -d "{}"'.format(';'.join(args.keys()), ';'.join(args.values()))
        print cmd
        os.system(cmd)


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test_denovo_ref_genome'),
                    TestFunction('test_denovo_ref_database'),
                    TestFunction('test_denovo_upload'),
                    TestFunction('test_upload_denovo')])
    unittest.TextTestRunner(verbosity=2).run(suite)

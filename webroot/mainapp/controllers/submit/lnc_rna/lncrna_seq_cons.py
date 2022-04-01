# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

from mainapp.controllers.project.lnc_rna_controller import LncRnaController
from mainapp.libs.signature import check_sig
import web
import datetime
import json
import unittest
import os

class LncrnaSeqConsAction(LncRnaController):
    '''
    last_modify: 2019.04.24
    '''
    def __init__(self):
        super(LncrnaSeqConsAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        args = ['task_id', 'submit_location', 'task_type', 'old_genome', 'new_genome', 'sequence']
        for arg in args:
            if not hasattr(data, arg):
                info = {'success': False, 'info': 'lack argument -> {}'.format(arg)}
                return json.dumps(info)
        if hasattr(data, 'lncrna_id') and data.sequence == 'single':
            lncrna_id = data.lncrna_id
            param_dict = {'lncrna_id': lncrna_id}
            options = {'keep_list': lncrna_id}
            to_files = ['lnc_rna.export_single_lncrna_list(keep_list)']
        elif hasattr(data, 'geneset_id') and data.sequence == 'multiple':
            geneset_id = data.geneset_id
            param_dict = {'geneset_id': geneset_id}
            options = {'keep_list': geneset_id}
            to_files = ['lnc_rna.export_lot_of_lncrna_list(keep_list)']

        task_id = data.task_id
        task_info = self.lnc_rna.get_task_info(task_id=task_id)
        project_sn = task_info['project_sn']
        time_now = datetime.datetime.now()
        created_ts = time_now.strftime('%Y-%m-%d %H:%M:%S')
        name = 'Lncrna_seq_cons_{}'.format(time_now.strftime('%Y%m%d_%H%M%S'))
        desc = 'lncRNA sequence conservation main table'
        param_dict.update({
            'task_id': task_id,
            'submit_location': data.submit_location,
            'task_type': int(data.task_type),
            'old_genome': data.old_genome,
            'new_genome': data.new_genome,
            'sequence': data.sequence
        })
        params = json.dumps(param_dict, sort_keys=True, separators=(',', ':'))
        main_info = {
            'task_id': task_id,
            'project_sn': project_sn,
            'created_ts': created_ts,
            'name': name,
            'desc': desc,
            'params': params,
            'status': 'start'
        }
        main_id = self.lnc_rna.insert_main_table('sg_lncrna_seq_cons', main_info)

        options.update({
            'lncrna_gtf': task_info['all_lnc_gtf'],
            'old_genome': data.old_genome,
            'new_genome': data.new_genome,
            'type_tsv': task_info['all_trans_type'],
            'main_id': str(main_id),
            'update_info': json.dumps({str(main_id): 'sg_lncrna_seq_cons'})
        })
        self.set_sheet_data(
            name='lnc_rna.report.lncrna_seq_cons',
            options=options,
            main_table_name=name,
            module_type='workflow',
            to_file=to_files,
            project_sn=project_sn,
            task_id=task_id
        )

        run_info = super(LncrnaSeqConsAction, self).POST()
        run_info['content'] = {'ids': {'id': str(main_id), 'name': name}}
        return json.dumps(run_info)

class TestFunction(unittest.TestCase):
    '''
    This is test for the controller. Just run this script to do test.
    '''
    def test_lncrna_id(self):
        cmd = 'python /mnt/ilustre/users/sanger-dev/biocluster/bin/webapitest.py post'
        cmd += ' -fr no'
        cmd += ' -c client03'
        cmd += ' -b http://192.168.12.101:9090'
        cmd += ' s/lnc_rna/lncrna_seq_cons'
        args = {
            'task_id': 'tsg_33915',
            'task_type': '2',
            'submit_location': 'lncrna_seq_cons',
            'lncrna_id': 'TCONS_00001138',
            'old_genome': 'hg38',
            'new_genome': 'mm10',
            'sequence': 'single'
        }
        cmd += ' -n "{}" -d "{}"'.format(';'.join(args.keys()), ';'.join(args.values()))
        print cmd
        os.system(cmd)

    def test_geneset_id(self):
        cmd = 'python /mnt/ilustre/users/sanger-dev/biocluster/bin/webapitest.py post'
        cmd += ' -fr no'
        cmd += ' -c client03'
        cmd += ' -b http://192.168.12.101:9090'
        cmd += ' s/lnc_rna/lncrna_seq_cons'
        args = {
            'task_id': 'tsg_33915',
            'task_type': '2',
            'submit_location': 'lncrna_seq_cons',
            'geneset_id': '5cbfc1108876b7e41f8b456d',
            'old_genome': 'hg38',
            'new_genome': 'mm10',
            'sequence': 'multiple'
        }
        cmd += ' -n "{}" -d "{}"'.format(';'.join(args.keys()), ';'.join(args.values()))
        print cmd
        os.system(cmd)

if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test_lncrna_id'), TestFunction('test_geneset_id')])
    unittest.TextTestRunner(verbosity=2).run(suite)

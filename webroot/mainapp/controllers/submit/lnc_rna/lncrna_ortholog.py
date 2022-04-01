# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

from mainapp.controllers.project.lnc_rna_controller import LncRnaController
from mainapp.libs.signature import check_sig
import web
import datetime
import json
from biocluster.config import Config
import unittest
import os

class LncrnaOrthologAction(LncRnaController):
    '''
    last_modify: 2019.05.05
    '''
    def __init__(self):
        super(LncrnaOrthologAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        args = ['task_id', 'submit_location', 'task_type', 'species', 'evalue', 'identity', 'sequence']
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
        name = 'Lncrna_ortholog_{}'.format(time_now.strftime('%Y%m%d_%H%M%S'))
        desc = 'lncRNA ortholog main table'
        param_dict.update({
            'task_id': task_id,
            'submit_location': data.submit_location,
            'task_type': int(data.task_type),
            'species': data.species,
            'evalue': data.evalue,
            'identity': data.identity,
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
        main_id = self.lnc_rna.insert_main_table('sg_lncrna_ortholog', main_info)

        db_col = Config().get_mongo_client(mtype='ref_rna_v2', dydb_forbid=True)[Config().get_mongo_dbname('ref_rna_v2', dydb_forbid=True)]['sg_genome_db']
        _id = max(record['_id'] for record in db_col.find({'name': data.species}) if 'lnc_dir' in record)
        options.update({
            'lncrna_fa': task_info['all_lnc_fa'],
            'target_fa': db_col.find_one({'_id': _id})['lnc_dir'],
            'evalue': float(data.evalue),
            'identity': float(data.identity),
            'type_tsv': task_info['all_trans_type'],
            'main_id': str(main_id),
            'update_info': json.dumps({str(main_id): 'sg_lncrna_ortholog'})
        })
        to_files.append('lnc_rna.export_target_fa(target_fa)')
        self.set_sheet_data(
            name='lnc_rna.report.lncrna_ortholog',
            options=options,
            main_table_name=name,
            module_type='workflow',
            to_file=to_files,
            project_sn=project_sn,
            task_id=task_id
        )

        run_info = super(LncrnaOrthologAction, self).POST()
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
        cmd += ' s/lnc_rna/lncrna_ortholog'
        args = {
            'task_id': 'tsg_33915',
            'submit_location': 'lncrna_ortholog',
            'task_type': '2',
            'lncrna_id': 'TCONS_00008759',
            'species': 'Arabidopsis_thaliana',
            'evalue': '1e-5',
            'identity': '50',
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
        cmd += ' s/lnc_rna/lncrna_ortholog'
        args = {
            'task_id': 'tsg_33915',
            'task_type': '2',
            'submit_location': 'lncrna_ortholog',
            'geneset_id': '5cbfc1108876b7e41f8b456d',
            'species': 'Arabidopsis_thaliana',
            'evalue': '1e-5',
            'identity': '50',
            'sequence': 'multiple'
        }
        cmd += ' -n "{}" -d "{}"'.format(';'.join(args.keys()), ';'.join(args.values()))
        print cmd
        os.system(cmd)

if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test_lncrna_id'), TestFunction('test_geneset_id')])
    unittest.TextTestRunner(verbosity=2).run(suite)

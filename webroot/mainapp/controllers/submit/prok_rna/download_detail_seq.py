# -*- coding: utf-8 -*-
# __author__ = 'fwy'

from mainapp.controllers.project.prok_rna_controller import ProkRNAController
from mainapp.libs.signature import check_sig
import web
import datetime
import json
import unittest
import os
from bson.objectid import ObjectId


class DownloadDetailSeqAction(ProkRNAController):
    '''
    last_modify: 20210603
    '''
    def __init__(self):
        super(DownloadDetailSeqAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        args = ['task_id', 'submit_location', 'task_type', 'extract_info', 'geneset_id']
        for arg in args:
            if not hasattr(data, arg):
                info = {'success': False, 'info': 'lack argument -> %s', "variables":[arg], "code" : "C3400102"}
                return json.dumps(info)

        task_id = data.task_id
        task_info = self.prok_rna.get_task_info(task_id=task_id)
        collection = self.db['sg_geneset']
        result = collection.find_one({"task_id": task_id, "main_id": ObjectId(data.geneset_id)})
        geneset_name = result["name"]
        # geneset_name = self.denovo_rna_v2.get_genesetname(task_id=task_id,geneset_id = data.geneset_id)
        project_sn = task_info['project_sn']
        time_now = datetime.datetime.now()
        created_ts = time_now.strftime('%Y-%m-%d %H:%M:%S')
        name = '{}_{}'.format(geneset_name,time_now.strftime('%Y%m%d_%H%M%S'))
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
            'desc': 'extract detail sequence main table',
            'created_ts': created_ts,
            'params': params,
            'status': 'start',
            'version': 'v1'
        }
        main_id = self.prok_rna.insert_main_table('sg_seq_extract', main_info)
        new_task_id = self.prok_rna.get_new_id(task_id)
        main_table_data = {'run_id': new_task_id}
        output_path = task_info["assemble_fa"]
        seq_detail = os.path.dirname(output_path)
        options = {
            'seq_detail_dir': seq_detail + '/',
            'extract_info': data.extract_info,
            'geneset_extract':data.geneset_id,
            'task_id': data.task_id,
            'main_id': str(main_id),
            'main_table_data': main_table_data,
            'update_info': json.dumps({str(main_id): 'sg_seq_extract'})
        }
        to_files = [
            'prok_rna.export_gene_list(geneset_extract)'
        ]
        self.set_sheet_data(
            name='prok_rna.report.download_detail_seq',
            options=options,
            main_table_name=name,
            module_type='workflow',
            to_file=to_files,
            project_sn=project_sn,
            task_id=task_id,
            new_task_id=new_task_id,
        )

        run_info = super(DownloadDetailSeqAction, self).POST()
        run_info['content'] = {'ids': {'id': str(main_id), 'name': name}}
        return json.dumps(run_info)


class TestFunction(unittest.TestCase):
    '''
    This is test for the controller. Just run this script to do test.
    '''
    def test(self):
        cmd = 'python /mnt/ilustre/users/sanger-dev/biocluster/bin/webapitest.py '
        cmd += 'post '
        cmd += "-fr no "
        cmd += '-dbversion 1 '
        cmd += '-c {} '.format("client03")
        cmd += ' -b http://bcl.tsg.com'
        cmd += ' s/prok_rna/download_detail_seq '
        args = {
            'task_id': 'gbfl_1pln6ksbdhrmsc3qrd1u4m',
            'submit_location': 'extract_seq_detail',
            'task_type': '2',
            'geneset_id':"60c22e30a4e1af2e4c4a4c1a",
            'extract_info': json.dumps({"seq_type":["cds", "protein"]}).replace('"', '\\"'),
        }
        cmd += ' -n "{}" -d "{}"'.format(';'.join(args.keys()), ';'.join(args.values()))
        print cmd
        os.system(cmd)


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)


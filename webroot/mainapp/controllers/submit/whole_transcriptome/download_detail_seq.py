# -*- coding: utf-8 -*-
# __author__ = 'fwy'

from mainapp.controllers.project.whole_transcriptome_controller import WholeTranscriptomeController
from mainapp.libs.signature import check_sig
import web
import datetime
import json
import unittest
import os
from bson.objectid import ObjectId


class DownloadDetailSeqAction(WholeTranscriptomeController):
    '''
    last_modify: 20210603
    '''
    def __init__(self):
        super(DownloadDetailSeqAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        print(data)
        args = ['task_id', 'submit_location', 'task_type', 'extract_info','geneset_id','level', 'category']
        for arg in args:
            if not hasattr(data, arg):
                info = {'success': False, 'info': 'lack argument -> %s', "variables":[arg], "code" : "C3400102"}
                return json.dumps(info)

        task_id = data.task_id
        task_info = self.whole_transcriptome.get_task_info(task_id=task_id)
        collection = self.db['geneset']
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
        main_id = self.whole_transcriptome.insert_main_table('seq_extract', main_info)
        new_task_id = self.whole_transcriptome.get_new_id(task_id)
        main_table_data = {'run_id': new_task_id}
        output_dir = task_info["output"]
        seq_detail = os.path.join(output_dir, "other/gene_detail/seqdetail/")
        options = {
            'seq_detail_dir': seq_detail,
            'extract_info': data.extract_info,
            'geneset_extract':data.geneset_id,
            'level': data.level,
            'category': data.category,
            'task_id': data.task_id,
            'main_id': str(main_id),
            'main_table_data': main_table_data,
            'update_info': json.dumps({str(main_id): 'seq_extract'})
        }
        to_files = [
            'whole_transcriptome.geneset.export_gene_list(geneset_extract)'
        ]
        self.set_sheet_data(
            name='whole_transcriptome.report.download_detail_seq',
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
        cmd += ' s/whole_transcriptome/download_detail_seq '
        args = {
            'task_id': 'rirc_pcuj135br7lchp8ca59i4g',
            'submit_location': 'extract_seq_detail',
            'task_type': '2',
            'level':'G',
            'geneset_id':"60dbe2f0f6b9e43049d415a6",
            'extract_info': json.dumps({"seq_type":["cds", "pep"]}).replace('"', '\\"'),
            'category': 'mRNA'
        }
        cmd += ' -n "{}" -d "{}"'.format(';'.join(args.keys()), ';'.join(args.values()))
        print cmd
        os.system(cmd)


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)


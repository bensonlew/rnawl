# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

from mainapp.controllers.project.small_rna_controller import SmallRnaController
from mainapp.libs.signature import check_sig
import web
import datetime
import json
import unittest
import os
from bson.objectid import ObjectId
import random
import re

class GenesetSelfAction(SmallRnaController):
    def __init__(self):
        super(GenesetSelfAction, self).__init__(instant=True)

    @check_sig
    def POST(self):
        input_data = web.input()

        ret = self.check_params(input_data)
        if ret != True:
            return ret

        self.task_id = input_data.task_id
        self.project_sn = self.small_rna.get_task_info(self.task_id)['project_sn']
        self.time_now = datetime.datetime.now()

        main_id = self.create_main_table(input_data)
        if main_id:
            self.small_rna.update_db_record(table_name='sg_geneset', record_id=main_id, main_id=main_id)
        else:
            info = {'success': False, 'info': 'fail to create main table in sg_geneset'}
            return json.dumps(info)

        ret = self.prepare_workflow(input_data, str(main_id))
        if ret != True:
            info = {'success': False, 'info': 'fail to prepare workflow'}
            return json.dumps(info)

        task_info = super(GenesetSelfAction, self).POST()
        task_info['content'] = {'ids': {'id': str(main_id), 'name': input_data.name}}
        return json.dumps(task_info)

    def check_params(self, data):
        expected_args = ['task_id', 'submit_location' ,'task_type', 'trait_path', 'gene_type', 'name', 'file_id']
        for arg in expected_args:
            if not hasattr(data, arg):
                info = {'success': False, 'info': 'Lack argument: {}'.format(arg)}
                return json.dumps(info)
        if data['gene_type'].upper() not in ['M', 'G']:
            info = {'success': False, 'info': 'The value of gene_type must be M or G'}
            return json.dumps(info)
        m = re.match(r'^\w+://\S+/.+$', data['trait_path'])
        if not m:
            info = {'success': False, 'info': 'The value of trait_path must be s3 remote input path'}
            return json.dumps(info)
        return True

    def create_main_table(self, data):
        desc = 'Geneset_upload_at_{}'.format(self.time_now.strftime('%Y%m%d_%H%M%S'))
        created_ts = self.time_now.strftime('%Y-%m-%d %H:%M:%S')
        params_dict = dict([
            ('task_id', self.task_id),
            ('submit_location', data.submit_location),
            ('task_type', data.task_type),
            ('trait_path', data.trait_path),
            ('gene_type', data.gene_type),
            ('name', data.name),
            ('file_id', data.file_id),
        ])
        params = json.dumps(params_dict, sort_keys=True, separators=(',', ':'))
        main_info = dict([
            ('task_id', self.task_id),
            ('type', data.gene_type),
            ('name', data.name),
            ('desc', desc),
            # ('gene_length', gene_length), # woulb be updated at workflow
            ('is_use', 0),
            ('status', 'start'),
            ('project_sn', self.project_sn),
            ('created_ts', created_ts),
            ('params', params),
        ])
        main_id = self.small_rna.insert_main_table('sg_geneset', main_info)
        return main_id

    def prepare_workflow(self, data, geneset_id):
        name = 'small_rna.report.geneset_self'
        main_table_name = 'GenesetGoDag_{}'.format(self.time_now.strftime('%Y%m%d_%H%M%S'))
        task_id = self.task_id
        project_sn = self.project_sn
        module_type = 'workflow'
        to_file = ['small_rna.export_geneset_from_query(genes)']
        options = dict([
            ('genes', task_id),
            ('file_path', data.trait_path),
            ('name', data.name),
            ('gene_type', data.gene_type),
            ('geneset_id', geneset_id),
            ('update_info', json.dumps({geneset_id: 'sg_geneset'})),
        ])
        self.set_sheet_data(
            name=name,
            main_table_name=main_table_name,
            task_id=task_id,
            project_sn=project_sn,
            module_type=module_type,
            to_file=to_file,
            options=options,
        )
        return True

class TestFunction(unittest.TestCase):
    """
    This is test for the controllers. Just run this script to do test.
    """
    def test_this(self):
        cmd = 'python /mnt/ilustre/users/sanger-dev/biocluster/bin/webapitest.py '
        cmd += 'post '
        cmd += '-fr no '
        cmd += '-c {} '.format("client03")
        cmd += 'i/small_rna/geneset_self '
        cmd += '-b http://192.168.12.102:9090 '
        args = {
            'task_id': 'tsg_33087',
            'submit_location': 'geneset_upload',
            'task_type': '1',
            'trait_path': 's3://commonbucket/files/m_188/188_5c08b2758d630/tsg_33087/35d029620e9b05f285542801797ab880.txt',
            'gene_type': 'M',
            'name': 'upload_geneset',
            'file_id': str(random.randint(9000000, 10000000))
        }
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
        print cmd
        os.system(cmd)

if __name__ == '__main__':
    unittest.main()

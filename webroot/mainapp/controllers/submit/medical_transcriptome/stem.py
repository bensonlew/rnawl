# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

from mainapp.controllers.project.medical_transcriptome_controller import MedicalTranscriptomeController
from mainapp.libs.signature import check_sig
import web
import pandas as pd
from bson.objectid import ObjectId
import numpy as np
import datetime
import json
import unittest
import os
import urlparse
import re
from biocluster.config import Config


class StemAction(MedicalTranscriptomeController):
    '''
    last_modify: 2019.07.15
    '''
    def __init__(self):
        super(StemAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        project_type = 'medical_transcriptome'
        db = Config().get_mongo_client(mtype=project_type)[Config().get_mongo_dbname(project_type)]
        data = web.input()
        ## cluster这个参数在框架代码中使用，增加处理，避免冲突
        if data.method == 'K':
            data_re = web.data()
            data1 = urlparse.parse_qs(data_re)
            if 'sanger' in data.cluster:
                cluster = data1['cluster'][0]
            else:
                cluster = data.cluster
        args = ['task_id', 'submit_location', 'task_type', 'exp_id', 'level', 'pheno_file', 'method', 'pheno_id']
        for arg in args:
            if not hasattr(data, arg):
                info = {'success': False, 'info': 'lack argument (%s)', "variables":[arg], "code" : "C2903908"}
                return json.dumps(info)
        for arg in {'SCM': ['number', 'unit', 'significance'], 'K': ['cluster']}[data.method]:
            if not hasattr(data, arg):
                info = {'success': False, 'info': 'lack argument (%s)', "variables":[arg], "code" : "C2903909"}
                return json.dumps(info)

        if data.method == 'SCM' and hasattr(data, 'unit'):
            if float(data.unit) >= 5:
                info = {'success': False, 'info': 'argument "unit" should be no larger than 5 (%s)', "variables":[data.unit], "code" : "C2903910"}
                return json.dumps(info)

        inter_dir = self.create_tmp_dir(data.task_id, 'stem')
        pheno_file = self.download_from_s3(data.pheno_file, inter_dir=inter_dir)
        c = re.findall("[^\w\\t\\n\-\r]", open(pheno_file, 'r').read())
        if c:
            d = list(set(c))
            return json.dumps({'success': False, 'info': "表格中存在以下特殊字符：{}".format("".join(d))})
        df = pd.read_table(pheno_file)
        df = df.reindex([i for i in df.columns if 'Unnamed: ' not in i], axis=1)
        raw_columns = df.columns
        df = df.rename({i: i.capitalize() for i in df.columns}, axis=1)
        req_columns = ['Sample', 'Group', 'Order']
        if list(df.columns) != req_columns:
            info = {'success': False, 'info': 'current columns of pheno file is {}; require {}'.format(
                str(list(raw_columns)), str(req_columns))}
            return json.dumps(info)
        connect_exp = db['sg_exp']
        record_exp = connect_exp.find_one({'task_id': data.task_id, 'main_id': ObjectId(data.exp_id)})
        is_rmbe = str(record_exp['is_rmbe']).lower()
        if is_rmbe == 'false':
            exp_id = data.exp_id
        if is_rmbe == 'true':
            exp_id = str(record_exp['batch_main_id'])
        sp2exp = self.medical_transcriptome.db['sg_exp_detail'].find_one(
            {'exp_id': ObjectId(exp_id)}, {'_id': 0, 'exp_id': 0, 'is_new': 0, 'seq_id': 0}
        )
        ncs = [s for s in df['Sample'] if s not in sp2exp]
        if ncs:
            info = {'success': False, 'info': 'can not find sample (%s) from pheno file in exp matrix', "variables":[ncs], "code" : "C2903912"}
            return json.dumps(info)
        if df['Order'].dtype == int:
            arr = df['Order'].unique()
            if not all(np.arange(arr.min(), arr.max() + 1) == arr):
                info = {'success': False, 'info': 'column "Order" should contain continuous integer (%s)', "variables":[arr], "code" : "C2903913"}
                return json.dumps(info)
        else:
            info = {'success': False, 'info': 'column "Order" can contain only integer (%s)', "variables":[df['Order'].dtype], "code" : "C2903914"}
            return json.dumps(info)
        df.to_csv(pheno_file, sep='\t', index=False)
        connect = db['sg_specimen_group']
        record = connect.find_one({'task_id': data.task_id})
        samples = record['specimen_names']
        print samples
        category_names = record['category_names']
        group_dict = dict(zip(category_names, samples))
        group_dict = json.dumps(group_dict)
        task_id = data.task_id
        task_info = self.medical_transcriptome.get_task_info(task_id=task_id)
        project_sn = task_info['project_sn']
        time_now = datetime.datetime.now()
        created_ts = time_now.strftime('%Y-%m-%d %H:%M:%S')
        name = 'Stem_{}_{}'.format(data.method, time_now.strftime('%Y%m%d_%H%M%S'))
        desc = 'Stem main table'
        param_dict = {
            'task_id': task_id,
            'submit_location': data.submit_location,
            'task_type': int(data.task_type),
            'exp_id': data.exp_id,
            'level': data.level,
            'pheno_file': data.pheno_file,
            'method': data.method,
            'pheno_id': data.pheno_id
        }
        if hasattr(data, 'geneset_id'):
            param_dict.update({'geneset_id': data.geneset_id})
        for arg in {'SCM': ['number', 'unit', 'significance'], 'K': ['cluster']}[data.method]:
            param_dict.update({arg: getattr(data, arg)})
        params = json.dumps(param_dict, sort_keys=True, separators=(',', ':'))
        main_info = {
            'task_id': task_id,
            'project_sn': project_sn,
            'created_ts': created_ts,
            'name': name,
            'desc': desc,
            'params': params,
            'status': 'start',
            'level': data.level,
            'version': 'v1'
        }
        main_id = self.medical_transcriptome.insert_main_table('sg_stem', main_info)
        new_task_id = self.medical_transcriptome.get_new_id(data.task_id)
        main_table_data = {'run_id': new_task_id}

        options = {
            'matrix': exp_id + ';' + data.level + ';' + is_rmbe,
            'pheno': data.pheno_file,
            'geneset': data.geneset_id + "," + exp_id + "," + data.level,
            'method': data.method,
            'main_id': str(main_id),
            'main_table_data': main_table_data,
            'update_info': json.dumps({str(main_id): 'sg_stem'}),
            'group_dict': group_dict
        }
        if data.method == 'SCM':
            options.update({
                'number': int(data.number),
                'unit': int(data.unit),
                'significance': float(data.significance),
            })
        elif data.method == 'K':
            options.update({'clusters': int(cluster)})
        to_files = [
            'medical_transcriptome.export_masigpro_geneset_medical(geneset)',
            'medical_transcriptome.export_exp_matrix_new(matrix)'
        ]
        self.set_sheet_data(
            name='medical_transcriptome.report.stem',
            options=options,
            main_table_name=name,
            module_type='workflow',
            to_file=to_files,
            project_sn=project_sn,
            new_task_id=new_task_id,
            task_id=task_id
        )

        run_info = super(StemAction, self).POST()
        run_info['content'] = {'ids': {'id': str(main_id), 'name': name}}
        return json.dumps(run_info)

class TestFunction(unittest.TestCase):
    '''
    This is test for the controller. Just run this script to do test.
    '''
    def test_scm(self):
        cmd = 'python /mnt/ilustre/users/sanger-dev/biocluster/bin/webapitest.py post'
        cmd += ' -fr no'
        cmd += ' -c client03'
        cmd += ' -b http://bcl.tsg.com'
        cmd += ' s/medical_transcriptome/stem'
        args = {
            'level': 'G',
            'task_id': 'medical_transcriptome',
            'submit_location': 'stem',
            'task_type': '2',
            'exp_id': '5f50cacf17b2bf5a6c8bfd88',
            'geneset_id': 'All',
            'pheno_file': '/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/medical_transcriptome/stem/pheno.txt',
            'method': 'SCM',
            'number': '50',
            'unit': '2',
            'significance': '0.05',
            'pheno_id': '1234567',
        }
        cmd += ' -n "{}" -d "{}"'.format(';'.join(args.keys()), ';'.join(args.values()))
        print cmd
        os.system(cmd)

    def test_k(self):
        cmd = 'python /mnt/ilustre/users/sanger-dev/biocluster/bin/webapitest.py post'
        cmd += ' -fr no'
        cmd += ' -c client03'
        cmd += ' -b http://192.168.12.101:9090'
        cmd += ' s/ref_rna_v2/stem'
        args = {
            'level': 'G',
            'task_id': 'tsg_34423',
            'submit_location': 'stem',
            'task_type': '2',
            'exp_id': '5d0706f417b2bf455831147f',
            'geneset_id': 'All',
            'pheno_file': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/ref_rna_v2/stem/tsg_34423.pheno.txt',
            'method': 'K',
            'cluster': '10'
        }
        cmd += ' -n "{}" -d "{}"'.format(';'.join(args.keys()), ';'.join(args.values()))
        print cmd
        os.system(cmd)

if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test_scm')])
    unittest.TextTestRunner(verbosity=2).run(suite)

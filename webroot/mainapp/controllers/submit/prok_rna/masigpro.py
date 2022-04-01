# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

from mainapp.controllers.project.prok_rna_controller import ProkRNAController
from mainapp.libs.signature import check_sig
import web
import pandas as pd
from bson.objectid import ObjectId
import numpy as np
import datetime
import json
import unittest
import os,re
import urlparse
import re

class MasigproAction(ProkRNAController):
    '''
    last_modify: 2019.07.16
    '''
    def __init__(self):
        super(MasigproAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        ## cluster这个参数在框架代码中使用，增加处理，避免冲突
        data_re = web.data()
        data1 = urlparse.parse_qs(data_re)
        if 'sanger' in data.cluster:
            cluster = data1['cluster'][0]
        else:
            cluster = data.cluster
        args = ['task_id', 'submit_location', 'task_type', 'exp_id',
                'design', 'cluster', 'method', 'design_id', 'geneset_id']
        for arg in args:
            if not hasattr(data, arg):
                info = {'success': False, 'info': 'lack argument -> %s', "variables":[arg], "code" : "C2903706"}
                return json.dumps(info)

        inter_dir = self.create_tmp_dir(data.task_id, 'masigpro')
        design = self.download_from_s3(data.design, inter_dir=inter_dir)
        c = re.findall("[^\w\\t\\n\-\r]", open(design, 'r').read())
        if c:
            d = list(set(c))
            return json.dumps({'success': False, 'info': "表格中存在以下特殊字符：{}".format("".join(d))})

        df = pd.read_table(design)
        df = df.reindex([i for i in df.columns if 'Unnamed: ' not in i], axis=1)
        raw_columns = df.columns
        df = df.rename({i: i.capitalize() for i in df.columns}, axis=1)
        req_columns = ['Sample', 'Time', 'Replicate']
        if len(df.columns) < 4:
            info = {'success': False, 'info': 'columns number should greater than 3 (%s)', "variables":[len(df.columns)], "code" : "C2903707"}
            return json.dumps(info)
        if list(df.columns)[:3] != req_columns:
            info = {'success': False, 'info': 'current columns of pheno file is {}; require {}'.format(
                str(list(raw_columns)), str(req_columns))}
            return json.dumps(info)
        sp2exp = self.prok_rna.db['sg_exp_detail'].find_one(
            {'exp_id': ObjectId(data.exp_id)}, {'_id': 0, 'exp_id': 0, 'is_new': 0, 'seq_id': 0}
        )
        ncs = [s for s in df['Sample'] if s not in sp2exp]
        if ncs:
            # info = {'success': False, 'info': 'can not find sample (%s) from pheno file in exp matrix', "variables":[ncs[0]], "code" : "C2903709"}
            info = {'success': False, 'info': u'can not find sample {} from pheno file in exp matrix(时序样本信息表中的样本名不存在)'.format(";".join(ncs)),}
            return json.dumps(info)
        if df['Time'].dtype not in [int, float]:
            info = {'success': False, 'info': 'column "Time" should contain continuous number (%s)'.format(
                df['Time'].unique())}
            return json.dumps(info)
        if df['Replicate'].dtype == int:
            arr = df['Replicate'].unique()
            if not all(np.arange(arr.min(), arr.max() + 1) == arr):
                info = {'success': False,
                        'info': 'column "Replicate" should contain continuous integer ({})'.format(arr)}
                return json.dumps(info)
        else:
            info = {'success': False,
                    'info': 'column "Replicate" can contain only integer ({})'.format(df['Replicate'].dtype)}
        for i in list(df.columns)[3:]:
            if df[i].dtype != int:
                info = {'success': False,
                        'info': 'column "{}" can contain only integer ({})'.format(i,df[i].dtype)}
                return json.dumps(info)
        replicates1 = set(df[df[df.columns[3]] == 0][df.columns[2]])
        replicates2 = set(df[df[df.columns[3]] == 1][df.columns[2]])
        if replicates1 & replicates2:
            info = {'success': False,
                    'info': 'column "Replicate" can not duplicate in different groups, ({}).'.format(replicates1 & replicates2)}
            return json.dumps(info)
        df.to_csv(design, sep='\t', index=False)

        task_id = data.task_id
        task_info = self.prok_rna.get_task_info(task_id=task_id)
        project_sn = task_info['project_sn']
        time_now = datetime.datetime.now()
        created_ts = time_now.strftime('%Y-%m-%d %H:%M:%S')
        name = 'Masigpro_{}'.format(time_now.strftime('%Y%m%d_%H%M%S'))
        desc = 'Masigpro main table'
        param_dict = {
            'task_id': task_id,
            'submit_location': data.submit_location,
            'task_type': int(data.task_type),
            'exp_id': data.exp_id,
            # 'exp_level': data.exp_level,
            'design': data.design,
            'cluster': cluster,
            'method': data.method,
            'design_id': data.design_id,
            'geneset_id': data.geneset_id
        }
        params = json.dumps(param_dict, sort_keys=True, separators=(',', ':'))
        # exp_level = data.exp_level
        main_info = {
            'task_id': task_id,
            'project_sn': project_sn,
            'created_ts': created_ts,
            'version' :"v3.1",
            'name': name,
            'desc': desc,
            'params': params,
            # 'exp_level': exp_level,
            'cluster': range(1, int(cluster) + 1),
            'status': 'start'
        }
        main_id = self.prok_rna.insert_main_table('sg_masigpro', main_info)
        new_task_id = self.prok_rna.get_new_id(data.task_id)
        main_table_data = {'run_id': new_task_id}

        exp_type = self.prok_rna.db['sg_exp'].find_one({'main_id': ObjectId(data.exp_id)})['exp_type']
        options = {
            'matrix': data.exp_id,
            'geneset': data.geneset_id,
            'design': design,
            'cluster': cluster,
            'method': data.method,
            'exp_type': exp_type,
            'main_id': str(main_id),
            'main_table_data': main_table_data,
            'update_info': json.dumps({str(main_id): 'sg_masigpro'})
        }
        to_files = [
            'prok_rna.export_masigpro_geneset(geneset)',
            'prok_rna.export_masigpro_matrix(matrix)'
        ]
        self.set_sheet_data(
            name='prok_rna.report.masigpro',
            options=options,
            main_table_name=name,
            module_type='workflow',
            to_file=to_files,
            project_sn=project_sn,
            new_task_id=new_task_id,
            task_id=task_id
        )

        run_info = super(MasigproAction, self).POST()
        run_info['content'] = {'ids': {'id': str(main_id), 'name': name}}
        return json.dumps(run_info)

class TestFunction(unittest.TestCase):
    '''
    This is test for the controller. Just run this script to do test.
    '''
    def test_hclust(self):
        cmd = 'python /mnt/ilustre/users/sanger-dev/biocluster/bin/webapitest.py post'
        cmd += ' -fr no'
        cmd += ' -c client03'
        cmd += ' -b http://192.168.12.101:9090'
        cmd += ' s/prok_rna/masigpro'
        args = {
            'task_id': 'tsg_33555',
            'submit_location': 'masigpro',
            'task_type': '2',
            'exp_id': '5c944ebd17b2bf489dd8b5cf',
            'exp_level': 'T',
            'geneset_id': '5cdd2da291fc2197308b4567',
            'design': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/ref_rna_v2/masigpro/s3.design.txt',
            'cluster': '9',
            'method': 'hclust'
        }
        cmd += ' -n "{}" -d "{}"'.format(';'.join(args.keys()), ';'.join(args.values()))
        print cmd
        os.system(cmd)

    def test_kmeans(self):
        cmd = 'python /mnt/ilustre/users/sanger-dev/biocluster/bin/webapitest.py post'
        cmd += ' -fr no'
        cmd += ' -c client03'
        cmd += ' -b http://192.168.12.101:9090'
        cmd += ' s/prok_rna/masigpro'
        args = {
            'task_id': 'tsg_33555',
            'submit_location': 'masigpro',
            'task_type': '2',
            'exp_id': '5c944ebd17b2bf489dd8b5cf',
            'exp_level': 'T',
            'geneset_id': 'All',
            'design': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/ref_rna_v2/masigpro/s3.design.txt',
            'cluster': '9',
            'method': 'kmeans'
        }
        cmd += ' -n "{}" -d "{}"'.format(';'.join(args.keys()), ';'.join(args.values()))
        print cmd
        os.system(cmd)

    def test_mclust(self):
        cmd = 'python /mnt/ilustre/users/sanger-dev/biocluster/bin/webapitest.py post'
        cmd += ' -fr no'
        cmd += ' -c client03'
        cmd += ' -b http://192.168.12.101:9090'
        cmd += ' s/prok_rna/masigpro'
        args = {
            'task_id': 'tsg_33555',
            'submit_location': 'masigpro',
            'task_type': '2',
            'exp_id': '5c944ebd17b2bf489dd8b5cf',
            'exp_level': 'T',
            'geneset_id': '5cdd2da291fc2197308b4567',
            'design': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/ref_rna_v2/masigpro/s3.design.txt',
            'cluster': '7',
            'method': 'Mclust'
        }
        cmd += ' -n "{}" -d "{}"'.format(';'.join(args.keys()), ';'.join(args.values()))
        print cmd
        os.system(cmd)

if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test_mclust')])
    unittest.TextTestRunner(verbosity=2).run(suite)

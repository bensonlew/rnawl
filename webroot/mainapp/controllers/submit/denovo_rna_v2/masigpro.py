# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import datetime
import json
import os
import unittest
import urlparse

import numpy as np
import pandas as pd
import web
from bson.objectid import ObjectId

from mainapp.controllers.project.denovo_rna_v2_controller import DenovoRnaV2Controller
from mainapp.libs.signature import check_sig


class MasigproAction(DenovoRnaV2Controller):
    '''
    last_modify: 2019.08.13
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
        args = ['task_id', 'submit_location', 'task_type', 'exp_id', 'exp_level',
                'design', 'cluster', 'method', 'design_id', 'geneset_id']
        for arg in args:
            if not hasattr(data, arg):
                info = {'success': False, 'info': 'lack argument -> %s', "variables": [arg], "code": "C1603401"}
                return json.dumps(info)

        inter_dir = self.create_tmp_dir(data.task_id, 'masigpro')
        design = self.download_from_s3(data.design, inter_dir=inter_dir)

        df = pd.read_table(design)
        df = df.reindex([i for i in df.columns if 'Unnamed: ' not in i], axis=1)
        raw_columns = df.columns
        df = df.rename({i: i.capitalize() for i in df.columns}, axis=1)
        req_columns = ['Sample', 'Time', 'Replicate']
        if len(df.columns) < 4:
            info = {'success': False, 'info': 'columns number should greater than 3 (%s)',
                    "variables": [len(df.columns)], "code": "C2903707"}
            return json.dumps(info)
        if list(df.columns)[:3] != req_columns:
            info = {'success': False, 'info': 'current columns of pheno file is %s; require %s'.format(
                raw_columns, req_columns)}
            return json.dumps(info)
        sp2exp = self.denovo_rna_v2.db['sg_exp_detail'].find_one(
            {'exp_id': ObjectId(data.exp_id)}, {'_id': 0, 'exp_id': 0, 'is_new': 0, 'seq_id': 0}
        )
        ncs = [s for s in df['Sample'] if s not in sp2exp]
        if ncs:
            info = {'success': False, 'info': 'can not find sample (%s) from pheno file in exp matrix',
                    "variables": [ncs], "code": "C2903709"}
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
            return json.dumps(info)
        df.to_csv(design, sep='\t', index=False)

        task_id = data.task_id
        task_info = self.denovo_rna_v2.get_task_info(task_id=task_id)
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
            'exp_level': data.exp_level,
            'design': data.design,
            'cluster': cluster,
            'method': data.method,
            'design_id': data.design_id,
            'geneset_id': data.geneset_id
        }
        params = json.dumps(param_dict, sort_keys=True, separators=(',', ':'))
        exp_level = data.exp_level
        main_info = {
            'task_id': task_id,
            'project_sn': project_sn,
            'created_ts': created_ts,
            'name': name,
            'desc': desc,
            'params': params,
            'exp_level': exp_level,
            'cluster': range(1, int(cluster) + 1),
            'status': 'start'
        }
        main_id = self.denovo_rna_v2.insert_main_table('sg_masigpro', main_info)

        exp_type = self.denovo_rna_v2.db['sg_exp'].find_one({'main_id': ObjectId(data.exp_id)})['exp_type']
        options = {
            'matrix': data.exp_id,
            'geneset': data.geneset_id,
            'design': design,
            'cluster': cluster,
            'method': data.method,
            'exp_type': exp_type,
            'main_id': str(main_id),
            'update_info': json.dumps({str(main_id): 'sg_masigpro'})
        }
        to_files = [
            'denovo_rna_v2.export_masigpro_geneset(geneset)',
            'denovo_rna_v2.export_masigpro_matrix(matrix)'
        ]
        self.set_sheet_data(
            name='denovo_rna_v2.report.masigpro',
            options=options,
            main_table_name=name,
            module_type='workflow',
            to_file=to_files,
            project_sn=project_sn,
            task_id=task_id
        )

        run_info = super(MasigproAction, self).POST()
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

    def test_hclust(self):
        cmd = 'python /mnt/ilustre/users/sanger-dev/biocluster/bin/webapitest.py post'
        cmd += ' -fr no'
        cmd += ' -c client03'
        cmd += ' -b http://bcl.tsg.com'
        cmd += ' s/denovo_rna_v2/masigpro'
        args = {
            'task_id': 'tsg_33857',
            'submit_location': 'masigpro',
            'task_type': '2',
            'exp_id': '5cb64b7617b2bf6187ad5bfe',
            'exp_level': 'T',
            'geneset_id': 'All',
            'design': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/denovo_rna_v2/timeseries/tsg_33857.design.txt',
            'design_id': '987654321',
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
        cmd += ' -b http://bcl.tsg.com'
        cmd += ' s/denovo_rna_v2/masigpro'
        args = {
            'task_id': 'tsg_33857',
            'submit_location': 'masigpro',
            'task_type': '2',
            'exp_id': '5cb64bc717b2bf6187af159d',
            'exp_level': 'G',
            'geneset_id': '5d47b73217b2bf35cade6f7b',
            'design': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/denovo_rna_v2/timeseries/tsg_33857.design.txt',
            'design_id': '987654321',
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
        cmd += ' -b http://bcl.tsg.com'
        cmd += ' s/denovo_rna_v2/masigpro'
        args = {
            'task_id': 'tsg_33857',
            'submit_location': 'masigpro',
            'task_type': '2',
            'exp_id': '5cb64bc717b2bf6187af159d',
            'exp_level': 'G',
            'geneset_id': '5d47b73217b2bf35cade6f7b',
            'design': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/denovo_rna_v2/timeseries/tsg_33857.design.txt',
            'design_id': '987654321',
            'cluster': '9',
            'method': 'Mclust'
        }
        cmd += ' -n "{}" -d "{}"'.format(';'.join(args.keys()), ';'.join(args.values()))
        print cmd
        os.system(cmd)


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test_hclust'),
                    TestFunction('test_kmeans'),
                    TestFunction('test_mclust')])
    unittest.TextTestRunner(verbosity=2).run(suite)

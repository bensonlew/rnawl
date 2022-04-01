# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import datetime
import json
import os
import unittest
import urlparse
from bson.objectid import ObjectId
import numpy as np
import pandas as pd
import web
from bson.objectid import ObjectId

from mainapp.controllers.project.whole_transcriptome_controller import WholeTranscriptomeController
from mainapp.libs.signature import check_sig


class MasigproAction(WholeTranscriptomeController):
    '''
    last_modify: 2019.07.16
    '''

    def __init__(self):
        super(MasigproAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        print data
        ## cluster这个参数在框架代码中使用，增加处理，避免冲突
        data_re = web.data()
        data1 = urlparse.parse_qs(data_re)
        if 'sanger' in data.cluster:
            cluster = data1['cluster'][0]
        else:
            cluster = data.cluster
        args = ['task_id', 'submit_location', 'task_type', 'level',
                'design', 'cluster', 'method', 'design_id', 'geneset_id', 'category', 'exp_id']
        for arg in args:
            if not hasattr(data, arg):
                info = {'success': False, 'info': 'lack argument -> %s', "variables": [arg], "code": "C2903706"}
                return json.dumps(info)

        inter_dir = self.create_tmp_dir(data.task_id, 'masigpro')
        if os.path.exists(data.design):
            design = data.design
        else:
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
        # exp_id = self.whole_transcriptome.db['exp'].find_one({'task_id': data.task_id, "level": data.level})["main_id"]
        exp_dict = self.whole_transcriptome.db['exp'].find_one({'task_id': data.task_id, "level": data.level, 'main_id': ObjectId(data.exp_id)})
        is_rmbe = str(exp_dict['is_rmbe'])
        if 'is_rmbe' not in exp_dict or is_rmbe == 'false':
            exp_id = str(exp_dict['main_id'])
        if is_rmbe == 'true':
            exp_id = str(exp_dict['batch_main_id'])

        sp2exp = self.whole_transcriptome.db['exp_detail'].find_one(
            {'exp_id': ObjectId(exp_id)}, {'_id': 0, 'exp_id': 0, 'is_new': 0, 'seq_id': 0}
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

        connect = self.db['specimen_group']
        record = connect.find_one({'task_id': data.task_id, 'library': 'long'})
        group_id = str(record['main_id'])
        samples = record['specimen_names']
        category_names = record['category_names']
        group_dict = dict(zip(category_names, samples))
        group_dict = json.dumps(group_dict, sort_keys=True, separators=(',', ':'))

        task_id = data.task_id
        task_info = self.whole_transcriptome.get_task_info(task_id=task_id)
        project_sn = task_info['project_sn']
        time_now = datetime.datetime.now()
        created_ts = time_now.strftime('%Y-%m-%d %H:%M:%S')
        name = 'Masigpro_{}'.format(time_now.strftime('%Y%m%d_%H%M%S'))
        desc = 'Masigpro main table'
        param_dict = {
            'task_id': task_id,
            'submit_location': data.submit_location,
            'task_type': int(data.task_type),
            'level': data.level,
            'design': data.design,
            'cluster': cluster,
            'method': data.method,
            'design_id': data.design_id,
            'geneset_id': data.geneset_id,
            'category': data.category,
            'exp_id': data.exp_id
        }
        params = json.dumps(param_dict, sort_keys=True, separators=(',', ':'))
        main_info = {
            'task_id': task_id,
            'project_sn': project_sn,
            'created_ts': created_ts,
            'name': name,
            'desc': desc,
            'params': params,
            'level': data.level,
            'category': data.category,
            'cluster': range(1, int(cluster) + 1),
            'status': 'start',
            'version': 'v1.1'
        }
        main_id = self.whole_transcriptome.insert_main_table('masigpro', main_info)

        way = self.whole_transcriptome.db['exp'].find_one({'main_id': ObjectId(data.exp_id)})['way']
        options = {
            'matrix': exp_id + ';' + data.level + ';' + is_rmbe,
            'geneset': data.geneset_id + "," + data.category + "," + data.level + "," + exp_id,
            'design': design,
            'cluster': cluster,
            'method': data.method,
            'level': data.level,
            'way': way,
            'group_dict': group_dict,
            'main_id': str(main_id),
            'update_info': json.dumps({str(main_id): 'masigpro'})
        }
        to_files = [
            'whole_transcriptome.whole_transcriptome.export_masigpro_geneset(geneset)',
            'whole_transcriptome_v1_1.whole_transcriptome.export_exp_matrix_new(matrix)'
        ]
        self.set_sheet_data(
            name='whole_transcriptome.report.masigpro',
            options=options,
            main_table_name=name,
            module_type='workflow',
            to_file=to_files,
            project_sn=project_sn,
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
        cmd += ' -b http://bcl.tsg.com'
        cmd += ' s/whole_transcriptome/masigpro'
        args = {
            'task_id': 'tsg_36088',
            'submit_location': 'masigpro',
            'task_type': '2',
            'level': 'G',
            'geneset_id': '5dbba38617b2bf6d13a3a818',
            'design': '/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/whole_transcriptome/masigpro/s3.design.txt',
            "design_id": "280234137",
            'category': 'mRNA',
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
        cmd += ' s/whole_transcriptome/masigpro'
        args = {
            'task_id': 'tsg_36088',
            'submit_location': 'masigpro',
            'task_type': '2',
            'level': 'G',
            'geneset_id': '5dbba38617b2bf6d13a3a818',
            'design': '/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/whole_transcriptome/masigpro/s3.design.txt',
            "design_id": "280234137",
            'category': 'mRNA',
            'cluster': '4',
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
        cmd += ' s/whole_transcriptome/masigpro'
        args = {
            'task_id': 'tsg_36088',
            'submit_location': 'masigpro',
            'task_type': '2',
            'level': 'G',
            'geneset_id': '5dbba38617b2bf6d13a3a818',
            'design': '/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/whole_transcriptome/masigpro/s3.design.txt',
            "design_id": "280234137",
            'category': 'mRNA',
            'cluster': '7',
            'method': 'Mclust'
        }
        cmd += ' -n "{}" -d "{}"'.format(';'.join(args.keys()), ';'.join(args.values()))
        print cmd
        os.system(cmd)


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test_kmeans')])
    unittest.TextTestRunner(verbosity=2).run(suite)

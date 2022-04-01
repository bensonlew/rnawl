# -*- coding: utf-8 -*-
# __author__ = 'zoujiaxun'

import datetime
import json
import os
import unittest
from bson.objectid import ObjectId
import web

from mainapp.controllers.project.whole_transcriptome_controller import WholeTranscriptomeController
from mainapp.libs.signature import check_sig


class BatchAction(WholeTranscriptomeController):
    def __init__(self):
        super(BatchAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        args = ['task_id', 'submit_location', 'task_type']
        args += ['library', 'level', 'exp_id', 'has_batch', 'batch_method']
        for arg in args:
            if not hasattr(data, arg):
                info = {'success': False, 'info': 'lack argument -> {}'.format(arg)}
                return json.dumps(info)

        task_id = data.task_id
        library = data.library

        exp_info = self.whole_transcriptome.get_exp_params_info(data.exp_id, data.task_id)
        project_sn = exp_info["project_sn"]
        params_json = {
            "exp_id": data.exp_id,
            'has_batch': data.has_batch,
            'level':data.level,
            'task_id': data.task_id,
            'submit_location': data.submit_location,
            'task_type': int(data.task_type),
            'batch_method': data.batch_method,
            'library': data.library
            # 'is_rmbe':'false'
        }
        if hasattr(data, "batch_matrix"):
            params_json.update({
                "batch_matrix": data.batch_matrix
            })
        if hasattr(data, "file_id"):
            params_json.update({
                'file_id': data.file_id
            })
        # batch json
        params_json = json.dumps(params_json, sort_keys=True, separators=(',', ':'))
        # exp json
        params = exp_info['params']
        params_dict = json.loads(params)
        if data.library == 'long':
            submit_location = 'exp_mrna'
        if data.library == 'small':
            submit_location = 'exp_mirna'
        if data.library == 'circle':
            submit_location = 'exp_circrna'
        params_dict.update({'submit_location': submit_location})
        params_dict_json = json.dumps(params_dict, sort_keys=True, separators=(',', ':'))
        connect = self.db['specimen_group']
        record = connect.find_one({'task_id': data.task_id, 'library': data.library})
        group_id = str(record['main_id'])
        samples = record['specimen_names']
        category_names = record['category_names']
        group_dict = dict(zip(category_names, samples))
        eps = 0
        for a in group_dict:
            if len(group_dict[a]) <3:
                eps += 1
            else:
                continue
        if eps == 0:
            ellipse = 'yes'
        else:
            ellipse = 'no'
        #判断批次表是否有误
        # if hasattr(data, "batch_matrix"):
        #     group_batch_dict = dict()
        #     for group in group_dict:
        #         for sample in group_dict[group]:
        #             if group not in group_batch_dict:
        #                 group_batch_dict[group] = [sample_batch_dict[sample]]
        #             else:
        #                 group_batch_dict[group].append(sample_batch_dict[sample])
        #     num = 0
        #     for group in group_batch_dict:
        #         if len(set(group_batch_dict[group])) == 1:
        #             num += 1
        #     if num == len(group_batch_dict):
        #         info = {'success': False, 'info': '批次表内容有误，所有组的组内不存在不同批次，无法做批次效应分析'}
        #         return json.dumps(info)
        group_dict = json.dumps(group_dict)
        connect_exp = self.db['exp']
        #record_exp = connect_exp.find_one({'task_id': data.task_id, 'main_id': ObjectId(data.exp_id), 'is_rmbe': False})
        record_exp = connect_exp.find_one({'task_id': data.task_id, 'main_id': ObjectId(data.exp_id)})
        try:
            batch = connect_exp.find_one({'task_id': data.task_id, 'main_id': ObjectId(data.exp_id), 'is_rmbe': 'true'})
            _id = batch['_id']
        except:
            pass
        way = record_exp['way']
        batch_version = record_exp['batch_version']
        connect_task = self.db['task']
        record_task = connect_task.find_one({'task_id': data.task_id})
        rnas = [rna for rna, lib in record_task['rna'].items() if lib == data.library]
        rnas_str = ";".join(rnas)
        name_batch = "ExpBatch" + '_' + data.level + '_'
        time_now = datetime.datetime.now()
        name_batch += time_now.strftime("%Y%m%d_%H%M%S")
        main_info_exp = dict(
            project_sn=project_sn,
            task_id=data.task_id,
            name=name_batch,
            level=data.level,
            way=way,
            version="v1.1",
            desc='{} exp_batch main table'.format(data.level),
            created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
            params=params_dict_json,
            status="start",
            batch_method=data.batch_method,
            group_dict=group_dict,
            is_rmbe='true',
            batch_version=batch_version,
            category=rnas_str,
            library=data.library
            # main_id=ObjectId(data.exp_id)
        )
        main_id_exp = self.whole_transcriptome.insert_main_table('exp', main_info_exp)

        name = 'Exp_batch_{}_{}'.format(library, time_now.strftime('%Y%m%d_%H%M%S'))
        task_info = self.whole_transcriptome.get_task_info(task_id=data.task_id)
        main_info_batch = {
            'task_id': data.task_id,
            'project_sn': task_info['project_sn'],
            'name': name,
            'created_ts': time_now.strftime('%Y-%m-%d %H:%M:%S'),
            'desc': 'exp batch main table',
            'params': params_json,
            'status': 'start',
            'library': library,
            'group_dict' : group_dict,
            'batch_version': batch_version,
            'exp_batch_id': ObjectId(str(main_id_exp)),
            'version': 'v1.1',
            'group_id': group_id,
            'level': data.level
        }
        if ellipse == 'yes':
            main_info_batch.update({
                'ellipse': 'yes'
            })
        if ellipse == 'no':
            main_info_batch.update({
                'ellipse': 'no'
            })
        main_id_batch = self.whole_transcriptome.insert_main_table('exp_batch', main_info_batch)

        sample_list = list()
        for sample in samples:
            sample_list += sample
        sample_str = ','.join(sample_list)
        options = {
            'has_batch': data.has_batch,
            # 'update_info': json.dumps({str(main_id_batch): "exp_batch", str(main_id_exp): 'exp'}),
            'update_info': json.dumps({str(main_id_batch): "exp_batch"}),
            # 'exp_id': str(main_id_exp),
            'exp_id': str(data.exp_id),
            'main_id_exp': str(main_id_exp),
            'sg_exp_batch_id': str(main_id_batch),
            'group_dict': group_dict,
            'batch_method': data.batch_method,
            'task_id': data.task_id,
            'params': params_json,
            'ellipse': ellipse,
            'exp_all': data.exp_id + ';' + sample_str + ';' + data.level + ';' + 'false',
            'level': data.level,
            'batch_version': batch_version,
            'library': library,
            'kind': 'all',
            'category': rnas_str,
            'is_rmbe': 'false'
        }
        if hasattr(data, "batch_matrix"):
            options.update({
                "batch_matrix": data.batch_matrix
            })

        to_files = ["whole_transcriptome.batch.export_exp_matrix_all(exp_all)"
                    ]


        self.set_sheet_data(
            name='whole_transcriptome_v1_1.report.batch',
            options=options,
            to_file=to_files,
            main_table_name=name,
            module_type='workflow',
            project_sn=project_sn,
            task_id=task_id
        )

        run_info = super(BatchAction, self).POST()
        run_info['content'] = {'ids': {'id': str(main_id_batch), 'name': name}}
        return json.dumps(run_info)
        if 'group_id' in data and str(data.group_id).lower() != 'all':
            _ = self.whole_transcriptome.update_group_is_use(data.task_id, data.group_id)
        if 'control_id' in data:
            _ = self.whole_transcriptome.update_group_compare_is_use(data.task_id, data.control_id)

        return json.dumps(run_info)


class TestFunction(unittest.TestCase):
    '''
    This is test for the controller. Just run this script to do test.
    '''

    def test(self):
        cmd = 'python /mnt/ilustre/users/sanger-dev/biocluster/bin/webapitest.py post'
        cmd += ' -fr no'
        cmd += ' -c client03'
        cmd += ' -b http://bcl.tsg.com'
        cmd += ' s/whole_transcriptome_v1_1/batch'
        args = {
            'task_id': 'tsg_38079',
            'submit_location': 'batch',
            'task_type': '2',
            'library': 'long',
            'level': 'G',
            'kind': 'all',
            'exp_id': "5f276f4417b2bf22947bf1f6",
            "has_batch": "True",
            "batch_method": 'combat',
            "batch_matrix": '/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/whole_transcriptome_v2/batch/batch_2'



        }
        cmd += ' -n "{}" -d "{}"'.format(';'.join(args.keys()), ';'.join(args.values()))
        print cmd
        os.system(cmd)


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)

# -*- coding: utf-8 -*-
# __author__ = 'konghualei,zoujiaxun'
from biocluster.config import Config
import datetime
import unittest
import json
import web
import os
from mainapp.controllers.project.ref_rna_v2_controller import RefRnaV2Controller
from mainapp.libs.signature import check_sig
from mainapp.models.mongo.ref_rna_v2 import *
from mbio.api.to_file.ref_rna_v2 import *
from bson.objectid import ObjectId

class BatchAction(RefRnaV2Controller):

    def __init__(self):
        super(BatchAction, self).__init__(instant=False)
    @check_sig
    def POST(self):
        project_type = 'ref_rna_v2'
        db = Config().get_mongo_client(mtype=project_type)[Config().get_mongo_dbname(project_type)]
        data = web.input()
        print data
        basic_args = ["task_id", "submit_location", 'task_type']
        basic_args += ["exp_id", 'exp_level', "has_batch", 'batch_method']

        # check arg
        for arg in basic_args:
            if not hasattr(data, arg):
                variables = list()
                variables.append(arg)
                info = {'success': False, 'info': "Lack argument: %s" % arg, 'code': 'C2901101', 'variables': variables}
                print info
                return json.dumps(info)

        inter_dir = self.create_tmp_dir(data.task_id, 'batch')
        if data.has_batch == "True":
            if not hasattr(data, "batch_matrix") or not data.batch_matrix:
                info = {'success': False, 'info': "请上传批次信息表"}
                return json.dumps(info)
        if hasattr(data, "batch_matrix"):
            batch_matrix = self.download_from_s3(data.batch_matrix, inter_dir=inter_dir)
            sample_batch_dict = dict()
            with open(batch_matrix, 'r') as f:
                for i in f.readlines():
                    if i.startswith('#'):
                        continue
                    else:
                        sample, batch = i.strip().split('\t')
                        # 注销部分是为了兼容一个报错项目，上传的group表中样本名以数字开头
                        # if sample == 'F4' and data.task_id == 'i-sanger_137095':
                        #     sample = '7F4'
                        sample_batch_dict[sample] = batch
            # 注销部分是为了兼容一个报错项目，上传的group表中样本名以数字开头
            # if data.task_id == 'i-sanger_137095':
            #     with open(batch_matrix, 'w') as new:
            #         new.write('#sample\tbatch\n')
            #         for i in sample_batch_dict:
            #             new.write(i+'\t'+sample_batch_dict[i]+'\n')

        exp_info = self.ref_rna_v2.get_exp_params_info(data.exp_id, data.task_id)
        project_sn = exp_info["project_sn"]
        params_json = {
            "exp_id": data.exp_id,
            'has_batch': data.has_batch,
            # 'exp_type': exp_info['exp_type'],
            'exp_level':data.exp_level,
            'task_id': data.task_id,
            'submit_location': data.submit_location,
            'task_type': int(data.task_type),
            'batch_method': data.batch_method

        }
        if hasattr(data, "batch_matrix"):
            params_json.update({
                "batch_matrix": data.batch_matrix
            })
        if hasattr(data, "batch_method"):
            params_json.update({
                "batch_method": data.batch_method
            })
        if hasattr(data, "file_id"):
            params_json.update({
                'file_id': data.file_id
            })

        params_json = json.dumps(params_json, sort_keys=True, separators=(',', ':'))
        params = exp_info['params']

        # if data.batch_method:
        #     name = "exp_batch" + '_' + data.batch_method + '_'
        # else:
        #     name = "exp_batch" + '_' + 'SVA' + '_'
        connect = db['sg_specimen_group']
        record = connect.find_one({'task_id': data.task_id})
        samples = record['specimen_names']
        category_names = record['category_names']
        group_dict = dict(zip(category_names, samples))


        eps = 0
        for a in group_dict:
            if len(group_dict[a]) <3:
                eps += 1
        # for sample_infos in group_dict.values():
        #     if len(sample_infos) < 3:
        #         eps += 1
            else:
                continue
        if eps == 0:
            ellipse = 'yes'
        else:
            ellipse = 'no'
        #判断批次表是否有误
        if hasattr(data, "batch_matrix"):
            group_batch_dict = dict()
            for group in group_dict:
                for sample in group_dict[group]:
                    if group not in group_batch_dict:
                        group_batch_dict[group] = [sample_batch_dict[sample]]
                    else:
                        group_batch_dict[group].append(sample_batch_dict[sample])
            num = 0
            for group in group_batch_dict:
                if len(set(group_batch_dict[group])) == 1:
                    num += 1
            if num == len(group_batch_dict):
                info = {'success': False, 'info': '批次表内容有误，所有组的组内不存在不同批次，无法做批次效应分析'}
                return json.dumps(info)

        group_dict = json.dumps(group_dict)

        connect_exp = db['sg_exp']
        record_exp = connect_exp.find_one({'task_id': data.task_id, 'main_id': ObjectId(data.exp_id)})
        print record_exp
        quant_method = record_exp['method']
        exp_type = record_exp['exp_type']
        libtype = record_exp['libtype']
        count_file = record_exp['count_file']
        name_batch = "ExpBatch" + '_' + data.exp_level + '_' + quant_method + '_' + exp_type + '_'
        time_now = datetime.datetime.now()
        name_batch += time_now.strftime("%Y%m%d_%H%M%S")
        main_info_exp = dict(
            project_sn=project_sn,
            task_id=data.task_id,
            name=name_batch,
            exp_level=data.exp_level,
            exp_type=exp_type.upper(),
            method=quant_method,
            version="v3.1",
            libtype=libtype,
            desc='{} exp_batch main table'.format(data.exp_level),
            created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
            params=params,
            status="start",
            count_file=count_file,
            batch_method=data.batch_method,
            group_dict=group_dict,
            raw_exp_id=ObjectId(data.exp_id)
        )
        main_id_exp = self.ref_rna_v2.insert_main_table('sg_exp', main_info_exp)

        name = 'exp_batch_{}'.format(time_now.strftime('%Y%m%d_%H%M%S'))
        task_info = self.ref_rna_v2.get_task_info(task_id=data.task_id)
        main_info_batch = {
            'task_id': data.task_id,
            'project_sn': task_info['project_sn'],
            'name': name,
            'created_ts': time_now.strftime('%Y-%m-%d %H:%M:%S'),
            'desc': 'exp batch main table',
            'params': params_json,
            'status': 'start',
            'group_dict': group_dict,
            # 'exp_batch_id': ObjectId(str(main_id_exp)),
            'relation_id': ObjectId(str(main_id_exp))
        }
        if ellipse == 'yes':
            main_info_batch.update({
                'ellipse': 'yes'
            })
        if ellipse == 'no':
            main_info_batch.update({
                'ellipse': 'no'
            })
        main_id_batch = self.ref_rna_v2.insert_main_table('sg_exp_batch', main_info_batch)
        new_task_id = self.ref_rna_v2.get_new_id(data.task_id)
        main_table_data = {'run_id': new_task_id}

        sample_list = list()
        for sample in samples:
            sample_list += sample
        sample_str = ','.join(sample_list)
        options = {
            "exp_matrix": data.exp_id + ';' + sample_str,
            'exp_other': data.exp_id + ';' + data.exp_level,
            "count": self.ref_rna_v2.get_count_path(data.exp_id),
            'has_batch': data.has_batch,
            'update_info': json.dumps({str(main_id_batch): "sg_exp_batch", str(main_id_exp): 'sg_exp'}),
            'exp_id': str(main_id_exp),
            'sg_exp_batch_id': str(main_id_batch),
            'group_table': data.task_id,
            # 'group_dict': group_dict
            'batch_method': data.batch_method,
            'task_id': data.task_id,
            'params':params_json,
            'ellipse': ellipse,
            'main_table_data': main_table_data,
        }
        if hasattr(data, "batch_matrix"):
            options.update({
                "batch_matrix": data.batch_matrix
            })
        # if hasattr(data, "batch_method"):
        #     options.update({
        #         "batch_method": data.batch_method
        #     })
        # prepare to file
        to_files = ["ref_rna_v2.export_exp_batch_matrix(exp_matrix)",
                    "ref_rna_v2.export_group_detail(group_table)",
                    "ref_rna_v2.export_exp_other_matrix(exp_other)"
                    ]
        task_name = 'ref_rna_v3.report.batch'

        self.set_sheet_data(name=task_name,
                            options=options,
                            main_table_name=name,  # 设置交互分析结果目录名
                            module_type="workflow",
                            to_file=to_files,
                            project_sn=project_sn,
                            new_task_id=new_task_id,
                            task_id=data.task_id)
        task_info = super(BatchAction, self).POST()
        task_info['content'] = {
            'ids': {
                'id': str(main_id_batch),
                'name': name
                }
        }
        # task_info['group_dict'] = group_dict
        return json.dumps(task_info)

        # 更新基因集的使用信息
        # self.whole_transcriptome.insert_geneset_info(data.geneset_id, 'geneset_cluster', str(main_id))
        if 'group_id' in data and str(data.group_id).lower() != 'all':
            _ = self.ref_rna_v2.update_group_is_use(data.task_id, data.group_id)
        if 'control_id' in data:
            _ = self.ref_rna_v2.update_group_compare_is_use(data.task_id, data.control_id)
        return json.dumps(task_info)


class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """

    def test_this(self):
        cmd = 'python /mnt/ilustre/users/sanger-dev/biocluster/bin/webapitest.py '
        cmd += 'post '
        cmd += "-fr no "
        cmd += '-c {} '.format("client03")
        cmd += "s/ref_rna_v3/batch "
        cmd += "-b http://bcl.tsg.com "
        args = dict(
            task_id="tsg_36994",

            task_type="2",
            submit_location="batch",
            exp_level='T',
            exp_id="5e5f039817b2bf64fcee6671",
            has_batch="True",
            batch_matrix='/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/batch_test/batch.txt',
            batch_method='combat'



        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)


if __name__ == '__main__':
    unittest.main()

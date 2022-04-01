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
        basic_args += ["exp_id", 'exp_level',"has_batch", 'batch_method']

        # check arg
        for arg in basic_args:
            if not hasattr(data, arg):
                variables = []
                variables.append(arg)
                info = {'success': False, 'info': "Lack argument: %s" % arg, 'code': 'C2901101', 'variables': variables}
                print info
                return json.dumps(info)

        exp_info = self.ref_rna_v2.get_exp_params_info(data.exp_id, data.task_id)
        project_sn = exp_info["project_sn"]
        params_json = {
            "exp_id": data.exp_id,
            'has_batch': data.has_batch,
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
        # if hasattr(data, "batch_method"):
        #     params_json.update({
        #         "batch_method": data.batch_method
        #     })

        params = json.dumps(params_json, sort_keys=True, separators=(',', ':'))


        # if data.batch_method:
        #     name = "exp_batch" + '_' + data.batch_method + '_'
        # else:
        #     name = "exp_batch" + '_' + 'SVA' + '_'
        connect = db['sg_specimen_group']
        record = connect.find_one({'task_id': data.task_id})
        samples = record['specimen_names']
        category_names = record['category_names']
        group_dict = dict(zip(category_names, samples))
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
        main_info = dict(
            project_sn=project_sn,
            task_id=data.task_id,
            name=name_batch,
            exp_level=data.exp_level,
            exp_type=exp_type.upper(),
            method=quant_method,
            version="v3",
            libtype=libtype,
            desc='{} exp_batch main table'.format(data.exp_level),
            created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
            params=params,
            status="start",
            count_file=count_file,
            batch_method = data.batch_method,
            group_dict=group_dict,
            raw_exp_id=ObjectId(data.exp_id)
        )





        # name = "exp_batch" + '_'
        # time_now = datetime.datetime.now()
        # name += time_now.strftime("%Y%m%d_%H%M%S")
        # main_info = dict(
        #     project_sn=project_sn,
        #     task_id=data.task_id,
        #     name=name,
        #     version="v1",
        #     exp_id=data.exp_id,
        #     created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
        #     desc='exp batch main table',
        #     params=params,
        #     status="start",
        #     group_dict=group_dict
        #
        # )
        # main_id = self.ref_rna_v2.insert_main_table('sg_exp_batch', main_info)
        main_id = self.ref_rna_v2.insert_main_table('sg_exp', main_info)

        # name_pca = "ExpPCA" + '_' + data.exp_level + '_'
        # time_now_pca = datetime.datetime.now()
        # name_pca += time_now_pca.strftime("%Y%m%d_%H%M%S")
        # main_info = dict(
        #     project_sn=project_sn,
        #     task_id=data.task_id,
        #     version="v3",
        #     name=name_pca,
        #     created_ts=time_now_pca.strftime('%Y-%m-%d %H:%M:%S'),
        #     exp_id=data.exp_id,
        #     desc='PCA main table',
        #     params=params,
        #     status="start"
        # )
        # main_id_pca = self.ref_rna_v2.insert_main_table('sg_exp_batch_pca', main_info)



        sample_list = list()
        for sample in samples:
            sample_list += sample
        sample_str = ','.join(sample_list)
        options = {
            "exp_matrix": data.exp_id + ';' + sample_str,
            'exp_other': data.exp_id + ';' + data.exp_level,
            "count": self.ref_rna_v2.get_count_path(data.exp_id),
            'has_batch': data.has_batch,
            'update_info': json.dumps({str(main_id): "sg_exp"}),
            'batch_id': str(main_id),
            'group_table': data.task_id,
            # 'group_dict': group_dict
            'batch_method': data.batch_method,
            'task_id': data.task_id

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
        task_name = 'tool_lab.batch'

        self.set_sheet_data(name=task_name,
                            options=options,
                            main_table_name=name_batch,  # 设置交互分析结果目录名
                            module_type="workflow",
                            to_file=to_files,
                            project_sn=project_sn,
                            task_id=data.task_id)
        task_info = super(BatchAction, self).POST()
        task_info['content'] = {
            'ids': {
                'id': str(main_id),
                'name': name_batch
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
        cmd += "s/tool_lab/batch "
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

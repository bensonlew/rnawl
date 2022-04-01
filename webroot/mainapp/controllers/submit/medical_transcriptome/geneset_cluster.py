#!/usr/bin/python
# -*- coding: utf-8 -*-
# __author__ = 'konghualei,20170426'

import web
import json
from mainapp.libs.signature import check_sig
import datetime
from mainapp.controllers.project.medical_transcriptome_controller import MedicalTranscriptomeController
from mbio.api.to_file.medical_transcriptome import *
from mainapp.models.mongo.medical_transcriptome import *
import os
import unittest
from collections import OrderedDict
from bson.objectid import ObjectId
from biocluster.config import Config


class GenesetClusterAction(MedicalTranscriptomeController):
    def __init__(self):
        super(GenesetClusterAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        basic_args = ["task_id", "submit_location", 'task_type']
        basic_args += ['group_id', 'group_dict', 'exp_id', 'level', 'use_group',
                       'geneset_id', 'sct', 'gct',
                       # 'scm', 'gcm',
                       # 'scd', 'gcd',
                       # 'n_clusters',
                       # 'geneset_type',
                       ]
        # check arg
        for arg in basic_args:
            if not hasattr(data, arg):
                variables = []
                variables.append(arg)
                info = {'success': False, 'info': "Lack argument: %s" % arg, 'code': 'C2901101', 'variables': variables}
                # print info
                return json.dumps(info)
        # if data.level != data.geneset_type:
        #     info = "level is not consistent with geneset_type!"
        #     return json.dumps(info)
        if data.group_id.lower() == 'all' and data.use_group.lower() == 'yes':
            info = {'success': False, 'info': "group 设置为all时,不能选择以分组均值分析", 'code': 'C2901102', 'variables': ''}
            return json.dumps(info)
        exp_info = self.medical_transcriptome.get_exp_params_info(data.exp_id, data.task_id)
        project_sn = exp_info["project_sn"]
        task_id = data.task_id
        group_dict = json.loads(data.group_dict, object_pairs_hook=OrderedDict)

        # create main table record
        level = exp_info['level']
        # exp_type, quant_method = exp_info['exp_type'], exp_info['method']

        params = dict(
            task_id=task_id,
            submit_location=data.submit_location,
            task_type=int(data.task_type),
            exp_id=data.exp_id,
            group_id=data.group_id,
            level=level,
            group_dict=group_dict,
            # quant_method=quant_method,
            use_group=data.use_group,
            geneset_id=data.geneset_id,
            # geneset_type=data.geneset_type,
            sct=data.sct,
            gct=data.gct,
        )

        n_clusters = 10
        scm = "complete"
        scd = "correlation"
        gcm = "average"
        gcd = "euclidean"
        if data.sct == 'no':
            pass
        elif data.sct == 'kmeans':
            params.update(dict(scd=data.scd, n_clusters=int(data.n_clusters), ))
            scd = data.scd
            n_clusters = int(data.n_clusters)
        else:
            params.update(dict(scm=data.scm, scd=data.scd))
            scm = data.scm
            scd = data.scd
        if data.gct.startswith('h'):
            params.update(dict(gcm=data.gcm, gcd=data.gcd, n_clusters=int(data.n_clusters)))
            gcm = data.gcm
            gcd = data.gcd
            n_clusters = int(data.n_clusters)
        elif data.gct.startswith('k'):
            params.update(dict(gcd=data.gcd, n_clusters=int(data.n_clusters)))
            gcd = data.gcd
            n_clusters = int(data.n_clusters)
        else:
            pass

        params = json.dumps(params, sort_keys=True, separators=(',', ':'))

        geneset_info = self.medical_transcriptome.get_main_info(data.geneset_id, 'sg_geneset', data.task_id)
        if not geneset_info:
            info = {"success": False, "info": "基因集不存在，请确认参数是否正确！!", 'variables': ''}
            return json.dumps(info)
        geneset_name = geneset_info['name']

        name = geneset_name + "_Cluster_"
        time_now = datetime.datetime.now()
        name += time_now.strftime("%Y%m%d_%H%M%S")
        main_info = dict(
            project_sn=project_sn,
            task_id=task_id,
            name=name,
            version="v1",
            group_dict=json.loads(data.group_dict, object_pairs_hook=OrderedDict),
            created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
            desc='geneset cluster main table',
            params=params,
            status="start",
            level=level,
        )
        main_id = self.medical_transcriptome.insert_main_table('sg_geneset_cluster', main_info)
        new_task_id = self.medical_transcriptome.get_new_id(task_id)
        main_table_data = {'run_id': new_task_id}

        # prepare option for workflow
        project_type = 'medical_transcriptome'
        db = Config().get_mongo_client(mtype=project_type)[Config().get_mongo_dbname(project_type)]
        connect_exp = db['sg_exp']
        record_exp = connect_exp.find_one({'task_id': data.task_id, 'main_id': ObjectId(data.exp_id)})
        if not record_exp:
            self.bind_object.set_error("意外错误，main_id:%s的背景基因在sg_exp中未找到！", variables=(data.exp_id,), )
        is_rmbe = str(record_exp['is_rmbe']).lower()
        if is_rmbe == 'false':
            exp_id = data.exp_id
        if is_rmbe == 'true':
            exp_id = str(record_exp['batch_main_id'])

        options = dict(
            exp_matrix=exp_id + ";" + data.geneset_id + ";" + data.level + ";" + is_rmbe,
            group_dict=data.group_dict,
            group=data.group_dict,
            cluster_main_id=str(main_id),
            n_clusters=n_clusters,
            group_id=data.group_id,
            use_group=data.use_group,
            scm=scm,
            scd=scd,
            sct=data.sct,
            gcm=gcm,
            gct=data.gct,
            gcd=gcd,
            gene_detail="|".join([task_id, data.level]),
            main_table_data=main_table_data,
            update_info=json.dumps({str(main_id): "sg_geneset_cluster"})  # to update sg_status
        )

        # prepare to file
        to_files = ["medical_transcriptome.export_geneset_exp_matrix_new(exp_matrix)",
                    "medical_transcriptome.export_group(group)",
                    "medical_transcriptome.get_gene_detail_new(gene_detail)"]

        # 把参数交给workflow运行相应的tool， 其中to_file用于准备tool的输入文件
        task_name = 'medical_transcriptome.report.geneset_cluster'
        self.set_sheet_data(name=task_name,
                            options=options,
                            main_table_name=name,  # 设置交互分析结果目录名
                            module_type="workflow",
                            to_file=to_files,
                            project_sn=project_sn,
                            new_task_id=new_task_id,
                            task_id=data.task_id)

        # 运行workflow 并传回参数
        task_info = super(GenesetClusterAction, self).POST()
        task_info['content'] = {
            'ids': {
                'id': str(main_id),
                'name': name
            }
        }
        # 更新基因集的使用信息
        self.medical_transcriptome.insert_geneset_info(data.geneset_id, 'sg_geneset_cluster', str(main_id))
        if 'group_id' in data and str(data.group_id).lower() != 'all':
            _ = self.medical_transcriptome.update_group_is_use(data.task_id, data.group_id)
        if 'control_id' in data:
            _ = self.medical_transcriptome.update_group_compare_is_use(data.task_id, data.control_id)
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
        cmd += "s/medical_transcriptome/geneset_cluster "
        cmd += "-b http://bcl.tsg.com "
        # cmd += "-b http://192.168.12.102:9090 "
        args = dict(
            task_id="medical_transcriptome",
            task_type="2",
            submit_location="genesetcluster",
            group_id="5f46228c17b2bf20e4e269e1",
            group_dict=r'{"H2": ["H1581_4", "H1581_5", "H1581_6"], "H3": ["H1581_7", "H1581_8", "H1581_9"], '
                       r'"S1": ["SNU16_1", "SNU16_2", "SNU16_3"], "H1": ["H1581_1", "H1581_2", "H1581_3"], '
                       r'"S3": ["SNU16_7", "SNU16_8", "SNU16_9"], "S2": ["SNU16_4", "SNU16_5", "SNU16_6"]}'.replace('"', '\\"'),
            exp_id="5f50cacf17b2bf5a6c8bfd88",
            level="G",
            geneset_id="5f45c70417b2bf78d9c9c174",
            quant_method='RSEM',
            use_group='no',
            n_clusters='10',
            sct="hierarchy",
            scm="complete",
            scd="correlation",
            gct="hierarchy",
            gcm="average",
            gcd="euclidean",
        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)


if __name__ == '__main__':
    unittest.main()

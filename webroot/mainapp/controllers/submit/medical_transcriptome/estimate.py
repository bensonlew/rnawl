#!/usr/bin/python
# -*- coding: utf-8 -*-
# __author__ = 'zoujiaxun'

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


class EstimateAction(MedicalTranscriptomeController):
    def __init__(self):
        super(EstimateAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        basic_args = ["task_id", "submit_location", 'task_type']
        basic_args += ['exp_file', 'id', 'platform', 'species', 'sct'
                       # 'scm', 'gcm',
                       # 'scd', 'gcd',
                       # 'n_clusters',
                       ]
        # check arg
        for arg in basic_args:
            if not hasattr(data, arg):
                variables = []
                variables.append(arg)
                info = {'success': False, 'info': "Lack argument: %s" % arg, 'code': 'C2901101', 'variables': variables}
                # print info
                return json.dumps(info)

        task_id = data.task_id
        task_info = self.medical_transcriptome.get_task_info(task_id=task_id)
        project_sn = task_info['project_sn']

        # create main table record
        params = dict(
            task_id=data.task_id,
            submit_location=data.submit_location,
            task_type=int(data.task_type),
            exp_file=data.exp_file,
            sct=data.sct,
        )
        if data.sct == 'no':
            pass
        elif data.sct == 'kmeans':
            params.update(dict(scd=data.scd, n_clusters=int(data.n_clusters), ))
        else:
            params.update(dict(scm=data.scm, scd=data.scd))


        params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        name = "Estimate" + '_'
        time_now = datetime.datetime.now()
        name += time_now.strftime("%Y%m%d_%H%M%S")
        main_info = dict(
            project_sn=project_sn,
            task_id=task_id,
            name=name,
            version="v1",
            created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
            desc='geneset cluster main table',
            params=params,
            status="start",
        )
        main_id = self.medical_transcriptome.insert_main_table('sg_tool_estimate', main_info)
        new_task_id = self.medical_transcriptome.get_new_id(task_id)
        main_table_data = {'run_id': new_task_id}

        # prepare option for workflow

        options = dict(
            exp_file=data.exp_file,
            id=data.id,
            platform=data.platform,
            species=data.species,
            sct=data.sct,
            scm=data.scm,
            scd=data.scd,
            main_table_data=main_table_data,
            main_id=str(main_id),
            update_info=json.dumps({str(main_id): "sg_tool_estimate"})  # to update sg_status
        )

        # prepare to file
        # to_files = ["medical_transcriptome.export_geneset_exp_matrix_new(exp_matrix)",
        #             "medical_transcriptome.export_group(group)",
        #             "medical_transcriptome.get_gene_detail_new(gene_detail)"]

        # 把参数交给workflow运行相应的tool， 其中to_file用于准备tool的输入文件
        task_name = 'medical_transcriptome.report.estimate'
        self.set_sheet_data(name=task_name,
                            options=options,
                            main_table_name=name,  # 设置交互分析结果目录名
                            module_type="workflow",
                            to_file=None,
                            project_sn=project_sn,
                            new_task_id=new_task_id,
                            task_id=data.task_id)

        # 运行workflow 并传回参数
        task_info = super(EstimateAction, self).POST()
        task_info['content'] = {
            'ids': {
                'id': str(main_id),
                'name': name
            }
        }
        # 更新基因集的使用信息
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
        cmd += "s/medical_transcriptome/estimate "
        cmd += "-b http://bcl.tsg.com "
        # cmd += "-b http://192.168.12.102:9090 "
        args = dict(
            task_id="medical_transcriptome",
            task_type="2",
            submit_location="estimate",
            exp_file='/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/tool_lab/nomogram/test2medical/count_test.txt',
            id='GeneSymbol',
            platform='illumina',
            species='Homo_sapiens',
            sct="hierarchy",
            scm="complete",
            scd="correlation",

        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)


if __name__ == '__main__':
    unittest.main()

# -*- coding: utf-8 -*-
import web
import json
import datetime
from mainapp.controllers.project.prok_rna_controller import ProkRNAController
from mbio.api.to_file.prok_rna import *
from mainapp.controllers.core.basic import Basic
from mainapp.libs.signature import check_sig
import unittest
import os
import re


class CheckGenomeAction(ProkRNAController):
    def __init__(self):
        super(CheckGenomeAction, self).__init__(instant=True)

    @check_sig
    def POST(self):
        data = web.input()
        basic_args = ["task_id", "submit_location", 'task_type']
        basic_args += ['seq_id', 'updownstream', 'tool_type']
        # check
        for arg in ['genome', 'gtf', 'in_type']:
            if not hasattr(data, arg):
                info = {'success': False, 'info': '缺少参数：{}'.format(arg)}
                return json.dumps(info)
            elif arg == "in_type" and data[arg] not in ["gtf", "gff"]:
                info = {'success': False, 'info': '参数in_type仅支持gtf/gff'}
                return json.dumps(info)
        time_now = datetime.datetime.now()
        task_id = time_now.strftime("%Y%m%d_%H%M%S")
        project_sn = time_now.strftime("%Y%m%d_%H%M%S")

        # create main table record
        params = dict(
            task_id=task_id,
            submit_location="vip_customer_center",
            task_type=1,
            genome = data.genome,
            gtf = data.gtf,
            in_type = data.in_type
        )
        params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        name = "CheckGenome" + '_' + time_now.strftime("%Y%m%d_%H%M%S")
        main_info = dict(
            project_sn=project_sn,
            task_id=task_id,
            name=name,
            created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
            desc='开始进行格式检查',
            params=params,
            status="start"
        )
        main_id = self.prok_rna.insert_main_table('sg_check_genome', main_info)
        # prepare option for workflow
        options = {
            "task_id": task_id,
            'params': params,
            "main_id": str(main_id),
            "genome": data.genome,
            "gtf": data.gtf,
            "in_type": data.in_type,
            "submit_location": "vip_customer_center",
            "update_info": json.dumps({str(main_id): "sg_check_genome"})  # to update sg_status
        }

        # 把参数交给workflow运行相应的tool， 其中to_file用于准备tool的输入文件
        task_name = 'prok_rna.report.check_genome'
        self.set_sheet_data_1(name=task_name,
                            options=options,
                            module_type="workflow",
                            project_sn=project_sn,
                            task_id=task_id,)

        # 运行workflow 并传回参数
        task_info = super(CheckGenomeAction, self).POST()
        task_info['content'] = {
            'ids': {
                'id': str(main_id),
                'name': name
                }
        }
        # 返回信息
        return json.dumps(task_info)


class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test_this(self):
        cmd = 'python /mnt/lustre/users/sanger-dev/wpm2/sanger_bioinfo/bin/webapitest.py '
        cmd += 'post '
        cmd += "-fr no "
        cmd += '-c {} '.format("client03")
        cmd += '-dbversion {} '.format(1)
        cmd += "i/prok_rna/check_genome "
        cmd += "-b http://wpm2.sanger.com "
        args = dict(
            genome="/mnt/lustre/users/sanger-dev/wpm2/workspace/20210722/Prokrna_2s9o_kmffc8ojen00hhqvbci5op/Download/output/GCF_000009345.1_ASM934v1/GCF_000009345.1_ASM934v1_genomic.fna",
            gtf="/mnt/lustre/users/sanger-dev/wpm2/workspace/20210722/Prokrna_2s9o_kmffc8ojen00hhqvbci5op/Download/output/GCF_000009345.1_ASM934v1/GCF_000009345.1_ASM934v1_genomic.gtf",
            in_type='gtf'
        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)


if __name__ == '__main__':
    unittest.main()

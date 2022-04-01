# -*- coding: utf-8 -*-
import web
import json
import datetime
from mainapp.controllers.project.lnc_rna_controller import LncRnaController
from mbio.api.to_file.lnc_rna import *
from mainapp.libs.signature import check_sig
import unittest
import os


class QuantAction(LncRnaController):
    def __init__(self):
        super(QuantAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        basic_args = ["task_id", "submit_location", 'task_type']
        basic_args += ['method', 'exp_type']
        # check arg
        for arg in basic_args:
            if not hasattr(data, arg):
                variables = []
                variables.append(arg)
                info = {'success': False, 'info': "Lack argument:  %s"% arg, 'code': 'C2902001', 'variables': variables}
                return json.dumps(info)
        sg_task = self.lnc_rna.get_task_info(data.task_id)
        project_sn = sg_task["project_sn"]
        libtype = self.lnc_rna.get_libtype(data.task_id)  # 随机选择一张主表获得libtype
        read_len = self.lnc_rna.get_mean_read_len(data.task_id)

        # create transcript main table record
        exp_level, exp_type, quant_method = 'T', data.exp_type.upper(), data.method
        name = "Exp" + '_' + exp_level + '_' + quant_method + '_' + exp_type + '_'
        time_now = datetime.datetime.now()
        name += time_now.strftime("%Y%m%d_%H%M%S")
        params = dict(
            task_id=data.task_id,
            submit_location=data.submit_location,
            task_type=int(data.task_type),
            method=data.method,
            exp_type=data.exp_type,
        )
        params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        main_info = dict(
            project_sn=project_sn,
            task_id=data.task_id,
            name=name,
            exp_level=exp_level,
            exp_type=exp_type.upper(),
            method=quant_method,
            libtype=libtype,
            desc='{} exp main table'.format(exp_level),
            created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
            params=params,
            status="start"
        )
        transcript_main_name = name
        transcript_main_id = self.lnc_rna.insert_main_table('sg_exp', main_info)

        # create transcript main table record
        exp_level, exp_type, quant_method = 'G', data.exp_type.upper(), data.method
        name = "Exp" + '_' + exp_level + '_' + quant_method + '_' + exp_type + '_'
        time_now = datetime.datetime.now()
        name += time_now.strftime("%Y%m%d_%H%M%S")
        params = dict(
            task_id=data.task_id,
            submit_location=data.submit_location,
            task_type=int(data.task_type),
            method=data.method,
            exp_type=data.exp_type,
        )
        params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        main_info = dict(
            project_sn=project_sn,
            task_id=data.task_id,
            name=name,
            exp_level=exp_level,
            exp_type=exp_type.upper(),
            method=quant_method,
            libtype=libtype,
            desc='{} exp main table'.format(exp_level),
            created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
            params=params,
            status="start"
        )
        gene_main_name = name
        gene_main_id = self.lnc_rna.insert_main_table('sg_exp', main_info)

        # prepare option for tool
        options = {
            "submit_location": data.submit_location,
            "task_type": data.task_type,
            "raw_task_id": data.task_id,
            "method": data.method,
            "libtype": libtype,
            "read_len": read_len,
            #"t2g": self.use_s3(sg_task["assemble_t2g"]),
            #"transcriptome": self.use_s3(sg_task["assemble_fa"]),
            #"fastq": sg_task['fastq'],
            "transcript_main_id": str(transcript_main_id),
            "gene_main_id": str(gene_main_id),
            "exp_type": exp_type,
            "task_id": data.task_id,
            "update_info": json.dumps({str(gene_main_id): "sg_exp", str(transcript_main_id): "sg_exp"})  # to update sg_status
        }

        # 把参数交给workflow运行相应的tool， 其中to_file用于准备tool的输入文件
        task_name = 'lnc_rna.report.quant'
        self.set_sheet_data(name=task_name,
                            options=options,
                            main_table_name=data.method + '_expression_quant',
                            module_type="workflow",
                            to_file=None,
                            project_sn=project_sn,
                            task_id=data.task_id)

        task_info = super(QuantAction, self).POST()
        task_info['content'] = {
            'ids': {
                'id': str(gene_main_id),
                'name': gene_main_name
                }
        }
        if 'group_id' in data and str(data.group_id).lower() != 'all':
            _ = self.lnc_rna.update_group_is_use(data.task_id, data.group_id)
        if 'control_id' in data:
            _ = self.lnc_rna.update_group_compare_is_use(data.task_id, data.control_id)
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
        cmd += "s/lnc_rna/quant "
        cmd += "-b http://192.168.12.101:9090 "
        args = dict(
            task_id="tsg_33924",
            task_type="2",
            submit_location="exp_deatil",
            method='RSEM',
            exp_type='FPKM',
        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(str(x) for x in arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)


if __name__ == '__main__':
    unittest.main()

# -*- coding: utf-8 -*-
import web
import json
import datetime
from collections import OrderedDict
from mainapp.controllers.project.ref_rna_controller import RefRnaController
from mbio.api.to_file.ref_rna import *
from mainapp.libs.signature import check_sig
import unittest
import os


class WgcnaModuleAction(RefRnaController):
    def __init__(self):
        super(WgcnaModuleAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        basic_args = ["task_id", "submit_location", 'task_type']
        basic_args += ['exp_level', 'wgcna_prepare_id', 'mergeCutHeight',
                       'power', 'minModuleSize', 'networkType', 'minKMEtoStay']
        # check arg
        for arg in basic_args:
            if not hasattr(data, arg):
                info = {'success': False, 'info': "Lack argument: {}".format(arg)}
                return json.dumps(info)
        exp_info = self.ref_rna.get_main_info(data.wgcna_prepare_id, 'sg_wgcna_prepare')
        project_sn = exp_info["project_sn"]
        task_id = data.task_id
        if data.power == "" or data.power == "NaN":
            prepare_info = self.ref_rna.get_power(data.wgcna_prepare_id)
            data.power = str(prepare_info["power_estimate"])

        # create main table record
        params = dict(
            task_id=task_id,
            submit_location=data.submit_location,
            task_type=int(data.task_type),
            wgcna_prepare_id=data.wgcna_prepare_id,
            exp_level=data.exp_level,
            mergeCutHeight=data.mergeCutHeight,
            power=data.power,
            minModuleSize=data.minModuleSize,
            networkType=data.networkType,
            minKMEtoStay=data.minKMEtoStay
        )
        params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        name = "WgcnaModule" + '_' + data.exp_level[0].upper() + '_'
        time_now = datetime.datetime.now()
        name += time_now.strftime("%Y%m%d_%H%M%S")
        main_info = dict(
            project_sn=project_sn,
            task_id=task_id,
            name=name,
            created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
            wgcna_prepare_id=data.wgcna_prepare_id,
            desc='wgcna module identification analysis main table',
            type=data.exp_level,
            params=params,
            status="start"
        )
        main_id = self.ref_rna.insert_main_table('sg_wgcna_module', main_info)

        # prepare option for workflow
        options = {
            "exp_matrix": data.wgcna_prepare_id,
            "main_id": str(main_id),
            "mergeCutHeight": data.mergeCutHeight,
            "power": data.power,
            "exp_level": data.exp_level,
            "networkType": data.networkType,
            "minModuleSize": data.minModuleSize,
            "minKMEtoStay": data.minKMEtoStay,
            "update_info": json.dumps({str(main_id): "sg_wgcna_module"})  # to update sg_status
        }
        # prepare to file
        to_files = ["ref_rna.export_wgcna_exp_matrix(exp_matrix)"]

        # 把参数交给workflow运行相应的tool， 其中to_file用于准备tool的输入文件
        task_name = 'ref_rna.report.wgcna_module'
        self.set_sheet_data(name=task_name,
                            options=options,
                            main_table_name=name,  # 设置交互分析结果目录名
                            module_type="workflow",
                            to_file=to_files,
                            project_sn=project_sn,
                            task_id=data.task_id)

        # 运行workflow 并传回参数
        task_info = super(WgcnaModuleAction, self).POST()
        task_info['content'] = {
            'ids': {
                'id': str(main_id),
                'name': name
                }
        }
        # task_info['group_dict'] = group_dict
        return json.dumps(task_info)


class TestFunction(unittest.TestCase):
    """
    This is test for the workflow. Just run this script to do test.
    """
    def test_this(self):
        cmd = 'python /mnt/ilustre/users/sanger-dev/biocluster/bin/webapitest.py '
        cmd += 'post '
        cmd += "-fr no "
        cmd += '-c {} '.format("client03")
        cmd += "s/ref_rna/wgcna_module "
        cmd += "-b http://192.168.12.102:9090 "
        args = dict(
            task_id="tsg_28226",
            task_type="2",
            submit_location="wgcna_module",
            wgcna_prepare_id="5ab85f79a4e1af4c15099262",
            exp_level="gene",
            mergeCutHeight="0.25",
            power="7",
            minModuleSize="30",
            networkType="signed",
            minKMEtoStay="0.5",
        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)


if __name__ == '__main__':
    unittest.main()

# -*- coding: utf-8 -*-

import datetime
import json
import os
import types
import unittest
from mainapp.controllers.project.medical_transcriptome_controller import MedicalTranscriptomeController
import web
from bson.objectid import ObjectId
from mbio.api.to_file.medical_transcriptome import *
from mainapp.libs.signature import check_sig


class WgcnaModuleAction(MedicalTranscriptomeController):
    def __init__(self):
        super(WgcnaModuleAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        basic_args = ["task_id", "submit_location", 'task_type']
        basic_args += ['level', 'wgcna_prepare_id', 'mergeCutHeight',
                       'power', 'minModuleSize', 'networkType', 'minKMEtoStay']
        # check arg
        for arg in basic_args:
            if not hasattr(data, arg):
                variables = []
                variables.append(arg)
                info = {'success': False, 'info': "Lack argument: %s" % arg, 'code': 'C2903201', 'variables': variables}
                return json.dumps(info)
        prepare_info = self.medical_transcriptome.get_main_info(data.wgcna_prepare_id, 'sg_wgcna_prepare', data.task_id)
        project_sn = prepare_info["project_sn"]
        task_id = data.task_id
        if data.power == "" or data.power == "NaN":
            prepare_info = self.medical_transcriptome.get_power(data.wgcna_prepare_id)
            data.power = str(prepare_info["power_estimate"])

        # create main table record
        params = dict(
            task_id=task_id,
            submit_location=data.submit_location,
            task_type=int(data.task_type),
            wgcna_prepare_id=data.wgcna_prepare_id,
            level=data.level,
            mergeCutHeight=data.mergeCutHeight,
            power=data.power,
            minModuleSize=data.minModuleSize,
            networkType=data.networkType,
            minKMEtoStay=data.minKMEtoStay
        )
        params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        name = "Module" + '_' + data.level[0].upper() + '_'

        time_now = datetime.datetime.now()
        name += time_now.strftime("%Y%m%d_%H%M%S")
        if isinstance(data.wgcna_prepare_id, types.StringTypes):
            wgcna_prepare_id = ObjectId(data.wgcna_prepare_id)
        elif isinstance(data.wgcna_prepare_id, ObjectId):
            wgcna_prepare_id = data.wgcna_prepare_id
        else:
            raise Exception("wgcna_prepare_id参数必须为字符串或者ObjectId类型!")
        main_info = dict(
            project_sn=project_sn,
            task_id=task_id,
            name=name,
            version="v1",
            created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
            wgcna_prepare_id=wgcna_prepare_id,
            desc='wgcna module identification analysis main table',
            level=data.level,
            params=params,
            status="start"
        )
        main_id = self.medical_transcriptome.insert_main_table('sg_wgcna_module', main_info)
        new_task_id = self.medical_transcriptome.get_new_id(data.task_id)
        main_table_data = {'run_id': new_task_id}

        # prepare option for workflow
        options = {
            "exp_matrix": data.wgcna_prepare_id,
            "main_id": str(main_id),
            "mergeCutHeight": data.mergeCutHeight,
            "power": data.power,
            "level": data.level,
            "networkType": data.networkType,
            "minModuleSize": data.minModuleSize,
            "minKMEtoStay": data.minKMEtoStay,
            'main_table_data': main_table_data,
            "update_info": json.dumps({str(main_id): "sg_wgcna_module"})  # to update status
        }
        # prepare to file
        to_files = ["medical_transcriptome.export_wgcna_exp_matrix(exp_matrix)"]

        # 把参数交给workflow运行相应的tool， 其中to_file用于准备tool的输入文件
        task_name = 'medical_transcriptome.report.wgcna_module'
        self.set_sheet_data(name=task_name,
                            options=options,
                            main_table_name=name,  # 设置交互分析结果目录名
                            module_type="workflow",
                            to_file=to_files,
                            project_sn=project_sn,
                            new_task_id=new_task_id,
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
        if 'group_id' in data and str(data.group_id).lower() != 'all':
            _ = self.medical_transcriptome.update_group_is_use(data.task_id, data.group_id)
        if 'control_id' in data:
            _ = self.medical_transcriptome.update_group_compare_is_use(data.task_id, data.control_id)
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
        cmd += "s/medical_transcriptome/wgcna_module "
        cmd += "-b http://bcl.tsg.com "
        args = dict(
            task_id="medical_transcriptome",
            task_type="2",
            submit_location="wgcna_module",
            wgcna_prepare_id="5f5830a817b2bf1d221b2cfe",
            level="G",
            mergeCutHeight="0.25",
            power="",
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

# -*- coding: utf-8 -*-
import web
import json
import datetime
from collections import OrderedDict
from mainapp.controllers.project.denovo_rna_v2_controller import DenovoRnaV2Controller
from mbio.api.to_file.denovo_rna_v2 import *
from mainapp.libs.signature import check_sig
import unittest
import os


class WgcnaPrepareAction(DenovoRnaV2Controller):
    def __init__(self):
        super(WgcnaPrepareAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        basic_args = ["task_id", "submit_location", 'task_type']
        basic_args += ['group_id', 'group_dict', 'exp_id', 'me', 'geneset_id', 'exp_level', 'cv']
        # check arg
        for arg in basic_args:
            if not hasattr(data, arg):
                var = []
                var.append(arg)
                info = {'success': False, 'info': "Lack argument: %s" % (arg), "code": 'C1602401', "variables": var}
                return json.dumps(info)
        exp_info = self.denovo_rna_v2.get_main_info(data.exp_id, 'sg_exp', data.task_id)
        project_sn = exp_info["project_sn"]
        task_id = data.task_id
        group_dict = json.loads(data.group_dict, object_pairs_hook=OrderedDict)
        # create main table record
        # exp_info = json.loads(exp_info['params'])
        exp_level, exp_type, quant_method = data.exp_level[0].upper(), exp_info['exp_type'], exp_info['method']
        params = dict(
            task_id=task_id,
            submit_location=data.submit_location,
            task_type=int(data.task_type),
            exp_id=data.exp_id,
            group_id=data.group_id,
            exp_level=data.exp_level,
            group_dict=group_dict,
            me=data.me,
            cv=data.cv,
            geneset_id=data.geneset_id
        )
        params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        name = "WgcnaPrepare" + '_' + exp_level + '_' + quant_method + '_' + exp_type.upper() + '_'
        time_now = datetime.datetime.now()
        name += time_now.strftime("%Y%m%d_%H%M%S")
        if type(params) == dict:
            params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        main_info = dict(
            project_sn=project_sn,
            task_id=task_id,
            name=name,
            created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
            exp_id=data.exp_id,
            desc='wgcna pre-processing analysis main table',
            type=data.exp_level,
            params=params,
            status="start"
        )
        main_id = self.denovo_rna_v2.insert_main_table('sg_wgcna_prepare', main_info)
        if data.geneset_id.lower() not in ["all", "refall", "none"]:
            self.denovo_rna_v2.insert_geneset_info(data.geneset_id, 'sg_wgcna_prepare', str(main_id))

        # prepare option for workflow
        if str(data.group_id).lower() == 'all':
            samples = group_dict['all']
            group_dict = OrderedDict([(x, [x]) for x in samples])
        options = {
            "exp_matrix": data.exp_id+";"+data.geneset_id,
            "group_dict": json.dumps(group_dict),
            "main_id": str(main_id),
            "me": data.me,
            "cv": data.cv,
            "exp_level": data.exp_level,
            "group_id": data.group_id,
            "update_info": json.dumps({str(main_id): "sg_wgcna_prepare"})  # to update sg_status
        }
        # prepare to file
        to_files = ["denovo_rna_v2.export_geneset_exp_matrix2(exp_matrix)"]

        # 把参数交给workflow运行相应的tool， 其中to_file用于准备tool的输入文件
        task_name = 'denovo_rna_v2.report.wgcna_prepare'
        self.set_sheet_data(name=task_name,
                            options=options,
                            main_table_name=name,  # 设置交互分析结果目录名
                            module_type="workflow",
                            to_file=to_files,
                            project_sn=project_sn,
                            task_id=data.task_id)

        # 运行workflow 并传回参数
        task_info = super(WgcnaPrepareAction, self).POST()
        task_info['content'] = {
            'ids': {
                'id': str(main_id),
                'name': name
                }
        }
        if 'group_id' in data and str(data.group_id).lower() != 'all':
            _ = self.denovo_rna_v2.update_group_is_use(data.task_id, data.group_id)
        if 'control_id' in data:
            _ = self.denovo_rna_v2.update_group_compare_is_use(data.task_id, data.control_id)
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
        cmd += "s/denovo_rna_v2/wgcna_prepare "
        cmd += "-b http://192.168.12.102:9090 "
        args = dict(
            task_id="denovo_rna_v2",
            task_type="2",
            exp_id='5a4c42a5a4e1af3d933435ae',
            group_id='59f3975f28fb4f19526ebee8',
            exp_level='gene',
            group_dict=json.dumps({"A": ["S1", "S2", "S3"], "B": ["S4", "S5", "S6"], "C": ["S7", "S8", "S9"]}).replace(
                '"', '\\"'),
            submit_location="wgcna_prepare",
            geneset_id='all',
            me="1",
            cv="0.2",
        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)


if __name__ == '__main__':
    unittest.main()

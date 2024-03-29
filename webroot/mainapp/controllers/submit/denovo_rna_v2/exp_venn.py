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


class ExpVennAction(DenovoRnaV2Controller):
    def __init__(self):
        super(ExpVennAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        basic_args = ["task_id", "submit_location", 'task_type']
        basic_args += ['group_id', 'group_dict', 'exp_id', 'threshold']
        # check arg
        for arg in basic_args:
            if not hasattr(data, arg):
                var = []
                var.append(arg)
                info = {'success': False, 'info': "Lack argument: %s"%(arg), "code": 'C1600801', "variables": var}
                return json.dumps(info)
        exp_info = self.denovo_rna_v2.get_exp_params_info(data.exp_id, data.task_id)
        project_sn = exp_info["project_sn"]
        task_id = data.task_id
        group_dict = json.loads(data.group_dict, object_pairs_hook=OrderedDict)
        # create main table record
        exp_level, exp_type, quant_method = exp_info['exp_level'], exp_info['exp_type'], exp_info['method']
        params = dict(
            task_id=task_id,
            submit_location=data.submit_location,
            task_type=int(data.task_type),
            exp_id=data.exp_id,
            group_id=data.group_id,
            exp_level=exp_level,
            group_dict=group_dict,
            # quant_method=quant_method,
            threshold=data.threshold,
        )
        params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        name = "ExpVenn" + '_' + exp_level + '_' + quant_method + '_' + exp_type.upper() + '_'
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
            version="v2",
            desc='Expression venn analysis main table',
            params=params,
            status="start"
        )
        main_id = self.denovo_rna_v2.insert_main_table('sg_exp_venn', main_info)

        # prepare option for workflow
        if str(data.group_id).lower() == 'all':
            samples = group_dict['all']
            group_dict = OrderedDict([(x, [x]) for x in samples])
        options = {
            "exp_matrix": data.exp_id,
            "group_dict": json.dumps(group_dict),
            "venn_main_id": str(main_id),
            "group": json.dumps(group_dict),
            "threshold": data.threshold,
            "update_info": json.dumps({str(main_id): "sg_exp_venn"})  # to update sg_status
        }
        # prepare to file
        to_files = ["denovo_rna_v2.export_exp_matrix(exp_matrix)",
                    "denovo_rna_v2.export_group(group)"]

        # 把参数交给workflow运行相应的tool， 其中to_file用于准备tool的输入文件
        task_name = 'denovo_rna_v2.report.exp_venn'
        self.set_sheet_data(name=task_name,
                            options=options,
                            main_table_name=name,  # 设置交互分析结果目录名
                            module_type="workflow",
                            to_file=to_files,
                            project_sn=project_sn,
                            task_id=data.task_id)

        # 运行workflow 并传回参数
        task_info = super(ExpVennAction, self).POST()
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
    This is test for the tool. Just run this script to do test.
    """
    def test_this(self):
        cmd = 'python /mnt/ilustre/users/sanger-dev/biocluster/bin/webapitest.py '
        cmd += 'post '
        cmd += "-fr no "
        cmd += '-c {} '.format("client03")
        cmd += "s/denovo_rna_v2/exp_venn "
        cmd += "-b http://bcl.tsg.com "
        args = dict(
            task_id="denovo_rna_v2_upgrade",
            task_type="2",
            submit_location="expvenn",
            group_id="5d47f30b17b2bf08c619bf78",
            group_dict=r'{"A":["A1", "A2", "A3"],"B": [ "B1", "B2", "B3"], "C": [ "C1", "C2", "C3"]}'.replace('"', '\\"'),
            exp_id="5d41458917b2bf10c8b3ffb0",
            exp_level="T",
            # geneset_id="5a35da8ba4e1af38e6fb83d9",
            # geneset_type="G",
            # quant_method='RSEM',
            threshold='2',
        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)


if __name__ == '__main__':
    unittest.main()

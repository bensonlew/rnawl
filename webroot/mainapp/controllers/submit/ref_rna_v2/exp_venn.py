# -*- coding: utf-8 -*-
import web
import json
import datetime
from collections import OrderedDict
from mainapp.controllers.project.ref_rna_v2_controller import RefRnaV2Controller
from mbio.api.to_file.ref_rna_v2 import *
from mainapp.libs.signature import check_sig
import unittest
import os


class ExpVennAction(RefRnaV2Controller):
    def __init__(self):
        super(ExpVennAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        basic_args = ["task_id", "submit_location", 'task_type']
        basic_args += ['group_id', 'group_dict', 'exp_id', 'threshold']
        basic_args += ['type']
        # check arg
        for arg in basic_args:
            if not hasattr(data, arg):
                variables = []
                variables.append(arg)
                info = {'success': False, 'info': "Lack argument: %s" % arg, 'code': 'C2900801', 'variables': variables}
                return json.dumps(info)
        exp_info = self.ref_rna_v2.get_exp_params_info(data.exp_id, data.task_id)
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
            type=data.type,
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
            version="v3",
            name=name,
            created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
            exp_id=data.exp_id,
            desc='Expression venn analysis main table',
            params=params,
            status="start"
        )
        main_id = self.ref_rna_v2.insert_main_table('sg_exp_venn', main_info)
        new_task_id = self.ref_rna_v2.get_new_id(task_id)
        main_table_data = {'run_id': new_task_id}

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
            "type": data.type,
            'main_table_data': main_table_data,
            "update_info": json.dumps({str(main_id): "sg_exp_venn"})  # to update sg_status
        }
        # prepare to file
        to_files = ["ref_rna_v2.export_exp_matrix(exp_matrix)",
                    "ref_rna_v2.export_group(group)"]

        # 把参数交给workflow运行相应的tool， 其中to_file用于准备tool的输入文件
        task_name = 'ref_rna_v2.report.exp_venn'
        self.set_sheet_data(name=task_name,
                            options=options,
                            main_table_name=name,  # 设置交互分析结果目录名
                            module_type="workflow",
                            to_file=to_files,
                            project_sn=project_sn,
                            new_task_id=new_task_id,
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
            _ = self.ref_rna_v2.update_group_is_use(data.task_id, data.group_id)
        if 'control_id' in data:
            _ = self.ref_rna_v2.update_group_compare_is_use(data.task_id, data.control_id)
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
        cmd += "s/ref_rna_v2/exp_venn "
        cmd += "-b http://192.168.12.101:9090 "
        args = dict(
            task_id="RefrnaV2_7320",
            task_type="2",
            submit_location="expvenn",
            exp_id='5b07ba2ca4e1af26ed3bc56a',
            group_id='5afd30bba4e1af301da5f379',
            exp_level='gene',
            threshold='2.0',
            group_dict=r'{"A1":["A1_1", "A1_2", "A1_3"],"A2": [ "A2_1", "A2_2", "A2_3"], "B1": [ "B1_1", "B1_2", "B1_3"], "B2": [ "B2_1", "B2_2", "B2_3"]}'.replace('"', '\\"'),
            type='ref',
        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(str(x) for x in arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTest(TestFunction('test_this'))
    unittest.TextTestRunner(verbosity=2).run(suite)

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


class ExpDistributionAction(RefRnaController):
    def __init__(self):
        super(ExpDistributionAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        basic_args = ["task_id", "submit_location", 'task_type']
        basic_args += ['group_id', 'exp_id', 'exp_level', 'group_dict']
        # check arg
        for arg in basic_args:
            if not hasattr(data, arg):
                info = {'success': False, 'info': "Lack argument: {}".format(arg)}
                return json.dumps(info)
        # exp_info = self.denovo_rna_v2.get_exp_params_info(data.exp_id, data.task_id)
        exp_info = self.ref_rna.get_main_info(data.exp_id, 'sg_express')
        project_sn = exp_info["project_sn"]
        task_id = data.task_id
        group_dict = json.loads(data.group_dict, object_pairs_hook=OrderedDict)
        if str(data.group_id).lower() == 'all':
            samples = group_dict['all']
            group_dict = OrderedDict([(x, [x]) for x in samples])

        # create main table record
        exp_level = data.exp_level
        exp_type = exp_info['name'].split('_')[2].lower()
        quant_method = exp_info['name'].split('_')[1]
        if quant_method.lower() == 'feacount':
            if exp_level != "gene":
                return json.dumps({'success': False, 'info': "You can only choose Gene for Analysis"})
        params = dict(
            task_id=task_id,
            submit_location=data.submit_location,
            task_type=int(data.task_type),
            exp_id=data.exp_id,
            group_id=data.group_id,
            exp_level=exp_level,
            group_dict=group_dict,
            # quant_method=quant_method,
        )
        if type(params) == dict:
            params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        name = "ExpDistribution" + '_' + exp_level[0].upper() + '_' + quant_method + '_'
        time_now = datetime.datetime.now()
        name += time_now.strftime("%Y%m%d_%H%M%S")
        main_info = dict(
            project_sn=project_sn,
            task_id=task_id,
            name=name,
            created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
            desc='density and box and volin plot main table',
            params=params,
            status="start"
        )
        main_id = self.ref_rna.insert_main_table('sg_exp_graph', main_info)

        # prepare option for workflow
        # if str(data.group_id).lower() == 'all':
        #     samples = group_dict['all']
        #     group_dict = OrderedDict([(x, [x]) for x in samples])
        options = {
            "exp_matrix": data.exp_id,
            # "group_dict": json.dumps(group_dict),
            "type": data.exp_level,  # gene/ transcript
            "express_level": exp_type,
            "graph_main_id": str(main_id),
            "update_info": json.dumps({str(main_id): "sg_exp_graph"})  # to update sg_status
        }

        # prepare to file
        to_files = ["ref_rna.export_express_matrix_level(exp_matrix)",]
        if str(data.group_id).lower() != 'all':
            options.update({"group_id": data.group_id, "group_detail": data.group_dict, })
            to_files.append("ref_rna.export_group_table_by_detail(group_id)")
        else:
            options.update({"group_id": None, "group_detail": data.group_dict, })

        # 把参数交给workflow运行相应的tool， 其中to_file用于准备tool的输入文件
        task_name = 'ref_rna.report.exp_distribution'
        self.set_sheet_data(name=task_name,
                            options=options,
                            main_table_name=name,  # 设置交互分析结果目录名
                            module_type="workflow",
                            to_file=to_files,
                            project_sn=project_sn,
                            task_id=data.task_id)

        # 运行workflow 并传回参数
        task_info = super(ExpDistributionAction, self).POST()
        task_info['content'] = {
            'ids': {
                'id': str(main_id),
                'name': name
                }
        }
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
        cmd += "s/ref_rna/exp_distribution "
        cmd += "-b http://192.168.12.102:9090 "
        args = dict(
            task_id="demo_sj_20171102_092507248",
            task_type="2",
            submit_location="exp_distr",
            group_id="59fa73f6a4e1af611f23a6a2",
            group_dict= r'{"A":["59fa73f3a4e1af611f239187","59fa73f3a4e1af611f239188","59fa73f3a4e1af611f239189"],"B":["59fa73f3a4e1af611f239182","59fa73f3a4e1af611f239183"]}'.replace('"', '\\"'),
            exp_id="59fa73f6a4e1af611f23a6ad",
            exp_level="gene",
        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)


if __name__ == '__main__':
    unittest.main()

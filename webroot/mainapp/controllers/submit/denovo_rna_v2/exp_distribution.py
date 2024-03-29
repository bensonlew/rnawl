# -*- coding: utf-8 -*-
import web
import json
import datetime
from collections import OrderedDict
from mainapp.controllers.project.denovo_rna_v2_controller import DenovoRnaV2Controller
from mbio.api.to_file.denovo_rna_v2 import *
from mainapp.libs.signature import check_sig
import os
import unittest

class ExpDistributionAction(DenovoRnaV2Controller):
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
                var = []
                var.append(arg)
                info = {'success': False, 'info': "Lack argument: %s"%(arg), "code": 'C1600601', "variables": var}
                return json.dumps(info)
        exp_info = self.denovo_rna_v2.get_exp_params_info(data.exp_id, data.task_id)
        project_sn = exp_info["project_sn"]
        task_id = data.task_id
        group_dict = json.loads(data.group_dict, object_pairs_hook=OrderedDict)

        # create main table record
        exp_level = exp_info['exp_level']
        exp_type, quant_method = exp_info['exp_type'], exp_info['method']
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
        name = "ExpDistribution" + '_' + exp_level + '_' + quant_method + '_'
        time_now = datetime.datetime.now()
        name += time_now.strftime("%Y%m%d_%H%M%S")
        main_info = dict(
            project_sn=project_sn,
            task_id=task_id,
            name=name,
            version="v2",
            created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
            desc='density and box and volin plot main table',
            params=params,
            status="start"
        )
        main_id = self.denovo_rna_v2.insert_main_table('sg_exp_graph', main_info)

        # prepare option for workflow
        if str(data.group_id).lower() == 'all':
            samples = group_dict['all']
            group_dict = OrderedDict([(x, [x]) for x in samples])
        options = {
            "exp_matrix": data.exp_id,
            "group_dict": json.dumps(group_dict),
            "graph_main_id": str(main_id),
            "update_info": json.dumps({str(main_id): "sg_exp_graph"})  # to update sg_status
        }

        # prepare to file
        to_files = ["denovo_rna_v2.export_exp_matrix(exp_matrix)",]

        # 把参数交给workflow运行相应的tool， 其中to_file用于准备tool的输入文件
        task_name = 'denovo_rna_v2.report.exp_distribution'
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
        if 'group_id' in data and str(data.group_id).lower() != 'all':
            _ = self.denovo_rna_v2.update_group_is_use(data.task_id, data.group_id)
        if 'control_id' in data:
            _ = self.denovo_rna_v2.update_group_compare_is_use(data.task_id, data.control_id)
        return json.dumps(task_info)



if __name__ == '__main__':
    class TestFunction(unittest.TestCase):
        """
        This is test for the tool. Just run this script to do test.
        """

        def test_this(self):
            cmd = 'python /mnt/ilustre/users/sanger-dev/biocluster/bin/webapitest.py '
            cmd += 'post '
            cmd += "-fr no "
            cmd += '-c {} '.format("client03")
            cmd += "s/denovo_rna_v2/exp_distribution "
            cmd += "-b http://bcl.tsg.com "
            args = dict(
                submit_location="exp_distribution",
                #db_type="ref_rna_v2",
                #geneset_id="5cc694d217b2bf5da04d69c3",
                #geneset_type="G",
                #anno_type="kegg",
                # kegg_table="5b1924afa4e1af33064178b3",
                exp_id="5d41458917b2bf10c8b3ffb0",
                group_id="5d47f30b17b2bf08c619bf78",
                exp_level='T',
                task_id="denovo_rna_v2_upgrade",
                task_type="2",
                group_dict=r'{"A":["A1", "A2", "A3"],"B": [ "B1", "B2", "B3"], "C": [ "C1", "C2", "C3"]}'.replace('"', '\\"'),
                #Draw_in_groups="yes"
                #type="origin",
                #method="bh"
            )

            arg_names, arg_values = args.keys(), args.values()
            cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
            print(cmd)
            os.system(cmd)


    unittest.main()
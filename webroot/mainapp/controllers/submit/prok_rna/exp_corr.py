# -*- coding: utf-8 -*-
import web
import json
import datetime
from collections import OrderedDict
from mainapp.controllers.project.prok_rna_controller import ProkRNAController
from mbio.api.to_file.prok_rna import *
from mainapp.libs.signature import check_sig
import unittest
import os


class ExpCorrAction(ProkRNAController):
    def __init__(self):
        super(ExpCorrAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        basic_args = ["task_id", "submit_location", 'task_type']
        basic_args += ['group_id', 'group_dict', 'exp_id', 'exp_level',
                       'corr_method', 'scm', 'scd', 'Draw_in_groups', 'log_base']
        # check arg# v3.2 exp_level 更新为ref_gene, novel_gene, sRNA
        for arg in basic_args:
            if not hasattr(data, arg):
                info = {'success': False, 'info': "Lack argument: {}".format(arg)}
                return json.dumps(info)
        exp_info = self.prok_rna.get_exp_params_info(data.exp_id, data.task_id)
        project_sn = exp_info["project_sn"]
        task_id = data.task_id
        group_dict = json.loads(data.group_dict, object_pairs_hook=OrderedDict)
        # create main table record
        exp_level = data.exp_level
        exp_type, quant_method = exp_info['exp_type'], exp_info['method']
        params = dict(
            task_id=task_id,
            submit_location=data.submit_location,
            task_type=int(data.task_type),
            exp_id=data.exp_id,
            group_id=data.group_id,
            exp_level=exp_level,
            group_dict=group_dict,
            scm=data.scm,
            scd=data.scd,
            # quant_method=quant_method,
            corr_method=data.corr_method,
            # type=data.type,
            Draw_in_groups=data.Draw_in_groups,
            log_base=data.log_base
        )
    

        params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        name = "ExpCorr" + '_' + exp_level + '_' + quant_method + '_'
        time_now = datetime.datetime.now()
        name += time_now.strftime("%Y%m%d_%H%M%S")
        main_info = dict(
            project_sn=project_sn,
            task_id=task_id,
            name=name,
            created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
            desc='correlation analysis main table',
            params=params,
            exp_level=exp_level,
            status="start"
        )
        main_id = self.prok_rna.insert_main_table('sg_exp_corr', main_info)

        # prepare option for workflow
        options = {
            "exp_matrix": data.exp_id,
            "group_dict": json.dumps(group_dict),
            "corr_main_id": str(main_id),
            "scm": data.scm,
            "scd": data.scd,
            "exp_level": exp_level,
            # "type": data.type,
            "corr_method": data.corr_method,
            "Draw_in_groups": data.Draw_in_groups,
            "log_base": str(data.log_base),
            "update_info": json.dumps({str(main_id): "sg_exp_corr"})  # to update sg_status
        }

        # prepare to file
        to_files = ["prok_rna.export_exp_matrix_prok(exp_matrix)",]

        # 把参数交给workflow运行相应的tool， 其中to_file用于准备tool的输入文件
        task_name = 'prok_rna.report.exp_corr'
        self.set_sheet_data(name=task_name,
                            options=options,
                            main_table_name=name,  # 设置交互分析结果目录名
                            module_type="workflow",
                            to_file=to_files,
                            project_sn=project_sn,
                            task_id=data.task_id)

        # 运行workflow 并传回参数
        task_info = super(ExpCorrAction, self).POST()
        task_info['content'] = {
            'ids': {
                'id': str(main_id),
                'name': name
                }
        }
        # task_info['group_dict'] = group_dict
        if 'group_id' in data and str(data.group_id).lower() != 'all':
            _ = self.prok_rna.update_group_is_use(data.task_id, data.group_id)
        if 'control_id' in data:
            _ = self.prok_rna.update_group_compare_is_use(data.task_id, data.control_id)
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
        cmd += "s/prok_rna/exp_corr "
        cmd += "-b http://192.168.12.102:9090 "
        args = dict(
            task_id="prok_rna_srna",
            task_type="2",
            submit_location="expcorr",
            exp_id="5b877be2a4e1af3eb5578ab9",
            group_id="5b7269fe77b3f3b1133b4a71",
            exp_level='mRNA',
            group_dict=r'{"WT": ["WT_1", "WT_2", "WT_3"], "rcsBKO": ["rcsBKO_1", "rcsBKO_2", "rcsBKO_3"]}'.replace('"', '\\"'),
            # type='ref',
            corr_method='pearson',
            scm='average',
            scd='euclidean',
        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(str(x) for x in arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTest(TestFunction('test_this'))
    unittest.TextTestRunner(verbosity=2).run(suite)

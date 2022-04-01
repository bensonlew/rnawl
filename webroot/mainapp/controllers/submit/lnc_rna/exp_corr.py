# -*- coding: utf-8 -*-
import web
import json
import datetime
from collections import OrderedDict
from mainapp.controllers.project.lnc_rna_controller import LncRnaController
from mbio.api.to_file.lnc_rna import *
from mainapp.libs.signature import check_sig
import unittest
import os


class ExpCorrAction(LncRnaController):
    def __init__(self):
        super(ExpCorrAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        basic_args = ["task_id", "submit_location", 'task_type','rna_type']#新增rna_type
        basic_args += ['group_id', 'group_dict', 'exp_id', 'exp_level',
                       'corr_method', 'scm', 'scd', 'type']
        # check arg
        for arg in basic_args:
            if not hasattr(data, arg):
                variables = []
                variables.append(arg)
                info = {'success': False, 'info': "Lack argument:  %s" % arg, 'code': 'C2900401', 'variables': variables}
                return json.dumps(info)
        exp_info = self.lnc_rna.get_exp_params_info(data.exp_id, data.task_id)
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
            scm=data.scm,
            scd=data.scd,
            # quant_method=quant_method,
            corr_method=data.corr_method,
            type=data.type,
            rna_type=data.rna_type
        )
        if hasattr(data, "log_base"):
            params.update({
                "log_base": data.log_base
            })

        params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        name = "ExpCorr" + '_' + exp_level + '_' + quant_method + '_'+data.rna_type+'_'
        time_now = datetime.datetime.now()
        name += time_now.strftime("%Y%m%d_%H%M%S")
        main_info = dict(
            project_sn=project_sn,
            task_id=task_id,
            name=name,
            created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
            desc='correlation analysis main table',
            params=params,
            status="start"
        )
        main_id = self.lnc_rna.insert_main_table('sg_exp_corr', main_info)

        # prepare option for workflow
        options = {
            "exp_matrix": data.exp_id,
            "group_dict": json.dumps(group_dict),
            "corr_main_id": str(main_id),
            "scm": data.scm,
            "scd": data.scd,
            "type": data.type,
            "rna_type": data.rna_type,
            "corr_method": data.corr_method,
            "update_info": json.dumps({str(main_id): "sg_exp_corr"})  # to update sg_status
        }

        if hasattr(data, "log_base"):
            options.update({
                "log_base": data.log_base
            })

        # prepare to file
        to_files = ["lnc_rna.export_exp_matrix(exp_matrix)",]

        # 把参数交给workflow运行相应的tool， 其中to_file用于准备tool的输入文件
        task_name = 'lnc_rna.report.exp_corr'
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
            _ = self.lnc_rna.update_group_is_use(data.task_id, data.group_id)
        if 'control_id' in data:
            _ = self.lnc_rna.update_group_compare_is_use(data.task_id, data.control_id)
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
        cmd += "s/lnc_rna/exp_corr "
        cmd += "-b http://192.168.12.101:9090 "
        args = dict(
            task_id="lnc_rna",
            task_type="2",
            submit_location="expcorr",
            exp_id='5c9da1ba17b2bf206a55d350',
            group_id='5c90920e17b2bf5bda46c152',
            exp_level='trans',
            group_dict=r'{"Con":["Con1","Con2","Con3","Con4","Con5"],"Vit": ["Vit1","Vit2","Vit3","Vit4","Vit5"]}'.replace('"', '\\"'),
            type='ref',
            corr_method='pearson',
            scm='average',
            scd='euclidean',
            rna_type="lncRNA",
        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(str(x) for x in arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTest(TestFunction('test_this'))
    unittest.TextTestRunner(verbosity=2).run(suite)

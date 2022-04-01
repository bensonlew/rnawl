# -*- coding: utf-8 -*-
import web
import json
import datetime
from collections import OrderedDict
from mainapp.controllers.project.itraq_and_tmt_controller import ItraqTmtController
from mbio.api.to_file.itraq_tmt import *
from mainapp.libs.signature import check_sig
import unittest

class ExpCorrAction(ItraqTmtController):
    def __init__(self):
        super(ExpCorrAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        basic_args = ["task_id", "submit_location", 'task_type']
        basic_args += ['group_id', 'group_dict', 'corr_method', 'scm', 'scd']
        # check arg
        for arg in basic_args:
            if not hasattr(data, arg):
                info = {'success': False, 'info': "Lack argument: %s", "variables":[arg], "code" : "C1900401"}
                return json.dumps(info)
            if arg.lower() == "null":
                info = {'success': False, 'info': "%s : is null or NULL", "variables":[arg], "code" : "C1900402"}
                return json.dumps(info)
        exp_info = self.itraq_tmt.get_exp_params_info_new(task_id=data.task_id, type="scaled")
        project_sn = exp_info["project_sn"]
        task_id = data.task_id
        group_dict = json.loads(data.group_dict, object_pairs_hook=OrderedDict)  # 需要保证样本顺序
        # create main table record
        params = dict(
            task_id=task_id,
            submit_location=data.submit_location,
            task_type=int(data.task_type),
            group_id=data.group_id,
            group_dict=group_dict,
            scm=data.scm,
            scd=data.scd,
            corr_method=data.corr_method,
        )
        params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        name = "SamCorr"
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
        main_id = self.itraq_tmt.insert_main_table('sg_express_corr', main_info)

        # prepare option for workflow
        options = {
            "exp_matrix": data.task_id,
            "group_dict": json.dumps(group_dict),
            "corr_main_id": str(main_id),
            "scm": data.scm,
            "scd": data.scd,
            "corr_method": data.corr_method,
            "update_info": json.dumps({str(main_id): "sg_express_corr"})  # to update sg_status
        }

        # prepare to file
        to_files = ["itraq_tmt.export_exp_matrix_scaled(exp_matrix)",]

        # 把参数交给workflow运行相应的tool， 其中to_file用于准备tool的输入文件
        task_name = 'itraq_and_tmt.report.exp_corr'
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
        cmd += "s/itraq_and_tmt/exp_corr "
        cmd += "-b http://192.168.12.102:9090 "
        args = dict(
            task_id="itraq_tmt",
            task_type="2",
            submit_location="expresscorr",
            group_id="59f3975f28fb4f19526ebee8",
            group_dict=json.dumps({"control_0":["control_0_1","control_0_2","control_0_3"],"BR_0":[ "BR_0_1","BR_0_2","BR_0_3" ], "control_3":["control_3_1","control_3_2","control_3_3"], "BR_3":["BR_3_1","BR_3_2","BR_3_3"]}).replace('"', '\\"'),
            scm="complete",
            scd="correlation",
            corr_method='pearson',
        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)


if __name__ == '__main__':
    unittest.main()

# -*- coding: utf-8 -*-
import web
import json
import datetime
from collections import OrderedDict
from mainapp.controllers.project.medical_transcriptome_controller import MedicalTranscriptomeController
from mbio.api.to_file.medical_transcriptome import *
from mainapp.libs.signature import check_sig
import unittest
import os


class ExpCorrAction(MedicalTranscriptomeController):
    def __init__(self):
        super(ExpCorrAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        basic_args = ["task_id", "submit_location", 'task_type']
        basic_args += ['group_id', 'group_dict', 'exp_id', 'level',
                       'corr_method', 'scm', 'scd', 'kind','Draw_in_groups']#20190521新增Draw_in_gropus
        # check arg
        for arg in basic_args:
            if not hasattr(data, arg):
                variables = []
                variables.append(arg)
                info = {'success': False, 'info': "Lack argument:  %s" % arg, 'code': 'C2900401', 'variables': variables}
                return json.dumps(info)
        exp_info = self.medical_transcriptome.get_exp_params_info(data.exp_id, data.task_id)
        is_rmbe = str(exp_info["is_rmbe"]).lower()
        if is_rmbe == 'false':
            exp_id = data.exp_id
        if is_rmbe == 'true':
            exp_id = str(exp_info['batch_main_id'])
        project_sn = exp_info["project_sn"]
        task_id = data.task_id
        group_dict = json.loads(data.group_dict, object_pairs_hook=OrderedDict)
        # create main table record
        level = exp_info['level']
        exp_type, quant_method = exp_info['exp_type'], exp_info['method']
        params = dict(
            task_id=task_id,
            submit_location=data.submit_location,
            task_type=int(data.task_type),
            exp_id=data.exp_id,
            group_id=data.group_id,
            level=level,
            group_dict=group_dict,
            scm=data.scm,
            scd=data.scd,
            # quant_method=quant_method,
            corr_method=data.corr_method,
            kind=data.kind,
            Draw_in_groups=data.Draw_in_groups,
            # is_rmbe=is_rmbe,
        )
        if hasattr(data, "log_base"):
            if data.log_base == "10":
                params.update({
                    "log_base": data.log_base
                })

        params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        name = "ExpCorr" + '_' + level + '_' + quant_method + '_'
        time_now = datetime.datetime.now()
        name += time_now.strftime("%Y%m%d_%H%M%S")
        main_info = dict(
            project_sn=project_sn,
            task_id=task_id,
            version="v1",
            name=name,
            created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
            desc='correlation analysis main table',
            params=params,
            status="start"
        )
        main_id = self.medical_transcriptome.insert_main_table('sg_exp_corr', main_info)
        new_task_id = self.medical_transcriptome.get_new_id(data.task_id)
        main_table_data = {'run_id': new_task_id}

        # prepare option for workflow
        options = {
            "exp_matrix": "{};{};{}".format(exp_id,level,is_rmbe),
            "group_dict": json.dumps(group_dict),
            "corr_main_id": str(main_id),
            "scm": data.scm,
            "scd": data.scd,
            "kind": data.kind,
            "corr_method": data.corr_method,
            "Draw_in_groups":data.Draw_in_groups,
            'main_table_data': main_table_data,
            "update_info": json.dumps({str(main_id): "sg_exp_corr"})  # to update sg_status
        }

        if hasattr(data, "log_base"):
            if data.log_base == "10":
                options.update({
                    "log_base": data.log_base
                })

        # prepare to file
        to_files = ["medical_transcriptome.export_exp_matrix_new(exp_matrix)",]

        # 把参数交给workflow运行相应的tool， 其中to_file用于准备tool的输入文件
        task_name = 'medical_transcriptome.report.exp_corr'
        self.set_sheet_data(name=task_name,
                            options=options,
                            main_table_name=name,  # 设置交互分析结果目录名
                            module_type="workflow",
                            to_file=to_files,
                            project_sn=project_sn,
                            new_task_id=new_task_id,
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
            _ = self.medical_transcriptome.update_group_is_use(data.task_id, data.group_id)
        if 'control_id' in data:
            _ = self.medical_transcriptome.update_group_compare_is_use(data.task_id, data.control_id)
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
        cmd += "s/medical_transcriptome/exp_corr "
        cmd += "-b http://bcl.tsg.com "
        args = dict(
            task_id="30v6_mbv5ll7kpn4ngt3vb15ste",
            task_type="2",
            submit_location="expcorr",
            exp_id='6019e89b17b2bf63f7d8b701', #主表有batch_main_id,is_rmbe=true
            # exp_id='5f50cb8917b2bf5a6c8d02bb', #主表无batch_main_id,is_rmbe=flase
            group_id='6019e63217b2bf63f7cfd33e',
            level='G',
            group_dict=r'{"shF2_87": ["shF2_3", "shF2_2", "shF2_1"], "shA1_87": ["shA1_3", "shA1_2", "shA1_1"], "ZC01_87": ["ZC01_3", "ZC01_2", "ZC01_1"], "DMSO_87": ["DMSO_3", "DMSO_2", "DMSO_1"], "YC49_87": ["YC49_3", "YC49_2", "YC49_1"], "shNS_87": ["shNS_3", "shNS_2", "shNS_1"]}'.replace('"','\\"'),
            kind='new',
            corr_method='pearson',
            scm='average',
            scd='euclidean',
            Draw_in_groups="yes",
            log_base="10",
        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(str(x) for x in arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTest(TestFunction('test_this'))
    unittest.TextTestRunner(verbosity=2).run(suite)

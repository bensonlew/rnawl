# -*- coding: utf-8 -*-
import web
import json
import datetime
from collections import OrderedDict
from mainapp.controllers.project.itraq_and_tmt_controller import ItraqTmtController
from mbio.api.to_file.itraq_tmt import *
from mainapp.libs.signature import check_sig
import unittest
import os


class WgcnaPrepareAction(ItraqTmtController):
    def __init__(self):
        super(WgcnaPrepareAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        basic_args = ["task_id", "submit_location", 'task_type']
        basic_args += ['group_id', 'group_dict', 'exp_id', 'me', 'proteinset_id', 'cv']
        # check arg
        for arg in basic_args:
            if not hasattr(data, arg):
                variables = []
                variables.append(arg)
                info = {'success': False, 'info': "Lack argument: %s" % arg, 'code': 'C2903501', 'variables': variables}
                return json.dumps(info)
        exp_info = self.itraq_tmt.get_main_info(data.exp_id, 'sg_express', data.task_id)
        project_sn = exp_info["project_sn"]
        task_id = data.task_id
        group_dict = json.loads(data.group_dict, object_pairs_hook=OrderedDict)
        # create main table record
        exp_info = json.loads(exp_info['params'])
        params = dict(
            task_id=task_id,
            submit_location=data.submit_location,
            task_type=int(data.task_type),
            exp_id=data.exp_id,
            group_id=data.group_id,
            group_dict=group_dict,
            me=data.me,
            cv=data.cv,
            proteinset_id=data.proteinset_id
        )
        params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        name = "WgcnaPrepare" + '_'
        time_now = datetime.datetime.now()
        name += time_now.strftime("%Y%m%d_%H%M%S")
        if type(params) == dict:
            params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        main_info = dict(
            project_sn=project_sn,
            task_id=task_id,
            name=name,
            version="v2",
            created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
            exp_id=data.exp_id,
            desc='wgcna pre-processing analysis main table',
            params=params,
            status="start"
        )
        main_id = self.itraq_tmt.insert_main_table('sg_wgcna_prepare', main_info)
        if data.proteinset_id.lower() not in ["all", "none"]:
            self.itraq_tmt.insert_proteinset_info(data.proteinset_id, 'sg_wgcna_prepare', str(main_id))

        # prepare option for workflow
        if str(data.group_id).lower() == 'all':
            samples = group_dict['all']
            group_dict = OrderedDict([(x, [x]) for x in samples])
        options = {
            "exp_matrix": data.exp_id+";"+data.proteinset_id,
            "group_dict": json.dumps(group_dict),
            "main_id": str(main_id),
            "me": data.me,
            "cv": data.cv,
            "group_id": data.group_id,
            "update_info": json.dumps({str(main_id): "sg_wgcna_prepare"})  # to update sg_status
        }
        # prepare to file
        to_files = ["itraq_tmt.export_proteinset_exp_matrix(exp_matrix)"]

        # 把参数交给workflow运行相应的tool， 其中to_file用于准备tool的输入文件
        task_name = 'itraq_and_tmt.report.wgcna_prepare'
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
        # task_info['group_dict'] = group_dict
        if 'group_id' in data and str(data.group_id).lower() != 'all':
            _ = self.itraq_tmt.update_group_is_use(data.task_id, data.group_id)
        if 'control_id' in data:
            _ = self.itraq_tmt.update_group_compare_is_use(data.task_id, data.control_id)
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
        cmd += "s/itraq_and_tmt/wgcna_prepare "
        cmd += "-b http://bcl.tsg.com "
        args = dict(
            task_id="tsg_34739",
            task_type="2",
            submit_location="wgcna_prepare",
            group_id="5d22a82717b2bf1c689db738",
            group_dict=r'{"F":["F_1","F_2","F_3"],"S":["S_1","S_2","S_3"]}'.replace(
                '"', '\\"'),
            exp_id="5d22a82817b2bf1c689dbed6",
            proteinset_id='all',
            me="0.1",
            cv="0.1",
        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)


if __name__ == '__main__':
    unittest.main()

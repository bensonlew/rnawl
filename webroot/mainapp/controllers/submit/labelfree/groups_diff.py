# -*- coding: utf-8 -*-
import web
import json
import datetime
from collections import OrderedDict
from mainapp.controllers.project.labelfree_controller import LabelfreeController
from mbio.api.to_file.labelfree import *
from mainapp.libs.signature import check_sig
import unittest
import os
from bson.objectid import ObjectId


class GroupsDiffAction(LabelfreeController):
    def __init__(self):
        super(GroupsDiffAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        basic_args = ["task_id", "submit_location", 'task_type']
        basic_args += ['group_id', 'group_dict', 'test_method', 'post_hoc', 'coverage']
        # type字段确认"origin"，"latest"
        # check arg
        for arg in basic_args:
            if not hasattr(data, arg):
                info = {'success': False, 'info': "Lack argument: %s", "variables":[arg], "code" : "C1900301"}
                return json.dumps(info)
            if arg.lower() == "null":
                info = {'success': False, 'info': "%s : is null or NULL", "variables":[arg], "code" : "C1900302"}
                return json.dumps(info)
        exp_info = self.labelfree.get_exp_params_info_new(task_id=data.task_id, type="ratio")
        project_sn = exp_info["project_sn"]
        task_id = str(data.task_id)


        group_dict = json.loads(data.group_dict)
        group_num = len(group_dict.keys())
        if group_num < 3:
            info = {'success': False, 'info': '分组数小于3组!'}
            return json.dumps(info)

        if not data.test_method in ["ow", 'kw']:  #ow : 单因素方差分析  kw： kw秩合检验
            info = {'success': False, 'info': '差异检验方法错误!'}
            return json.dumps(info)

        if not data.post_hoc in ["scheffe",'tukeykramer','gameshowell','welchuncorrected']:
          info = {'success': False, 'info': 'Post hoc 检验方法错误'}
          return json.dumps(info)


        # create main table record
        params = dict(
            task_id=task_id,
            submit_location=data.submit_location,
            task_type=int(data.task_type),
            group_id=data.group_id,
            group_dict=group_dict,
            test_method=data.test_method,
            post_hoc=data.post_hoc,
            coverage=data.coverage,
        )
        params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        name = "GroupsDiff" + '_'
        time_now = datetime.datetime.now()
        name += time_now.strftime("%Y%m%d_%H%M%S")
        if type(params) == dict:
            params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        main_info = dict(
            project_sn=project_sn,
            task_id=task_id,
            name=name,
            created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
            desc="多组比较和筛选",
            group_dict=group_dict,
            group_id=data.group_id,
            params=params,
            status="start",
        )

        main_id = self.labelfree.insert_main_table('sg_groups_diff', main_info)

        # prepare option for workflow
        options = {
            "coverage": data.coverage,
            "group": json.dumps(group_dict),
            "group_dict": json.dumps(group_dict),
            "group_name": "group",
            "main_table_id": str(main_id),
            "protein_table": task_id,
            "post_hoc": data.post_hoc,
            "test_method": data.test_method,
            "update_info": json.dumps({str(main_id): "sg_groups_diff"}),
        }
        # prepare to file
        to_files = ["labelfree.export_exp_matrix2(protein_table)",
                    "labelfree.export_group(group)"]

        # 把参数交给workflow运行相应的tool， 其中to_file用于准备tool的输入文件
        task_name = 'labelfree.report.groups_diff'
        self.set_sheet_data(name=task_name,
                            options=options,
                            main_table_name=name,  # 设置交互分析结果目录名
                            module_type="workflow",
                            to_file=to_files,
                            project_sn=project_sn,
                            task_id=task_id)

        # 运行workflow 并传回参数
        task_info = super(GroupsDiffAction, self).POST()
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
        cmd = 'python /mnt/lustre/users/sanger-dev/wpm2/sanger_bioinfo/bin/webapitest.py '
        cmd += 'post '
        cmd += "-fr no "
        cmd += '-c {} '.format("client03")
        cmd += '-dbversion {} '.format(1)
        cmd += "s/itraq_and_tmt/groups_diff "
        cmd += "-b http://wpm2.sanger.com "
        args = dict(
            task_id="jssn_svfguokp3g35cbuue5v31s",
            submit_location="groups_diff",
            task_type="2",
            group_id="611db7a8174ea0e30d956202",
            group_dict=json.dumps({"DMSO": ["DMSO_1", "DMSO_2"], "CG3556": ["CG3556_1", "CG3556_2"],"CG3926": ["CG3926_1", "CG3926_2"],"CG4209": ["CG4209_1", "CG4209_2"],}).replace('"', '\\"'),
            test_method="kw",
            post_hoc="scheffe",
            coverage="0.95"
        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)


if __name__ == '__main__':
    unittest.main()

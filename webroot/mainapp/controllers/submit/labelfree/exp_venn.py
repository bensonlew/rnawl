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

class ExpVennAction(LabelfreeController):
    def __init__(self):
        super(ExpVennAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        print(dict(data))
        basic_args = ["task_id", "submit_location", 'task_type']
        basic_args += ['group_id', 'group_dict']
        # check arg
        for arg in basic_args:
            if not hasattr(data, arg):
                info = {'success': False, 'info': "Lack argument: {}".format(arg)}
                return json.dumps(info)
        if len(json.loads(data.group_dict)) == len(json.loads(data.group_dict).values()):
            if len(json.loads(data.group_dict)) < 2 or len(json.loads(
                    data.group_dict)) > 6:
                info = {'success': False, 'info': "you should select less "
                                                  "than 7 or more than 1 "
                                                  "sample"}
                return json.dumps(info)

        exp_info = self.labelfree.get_exp_params_info(data.task_id)
        project_sn = exp_info["project_sn"]
        task_id = data.task_id
        group_dict = json.loads(data.group_dict, object_pairs_hook=OrderedDict)
        # create main table record
        params = dict(
            task_id=task_id,
            submit_location=data.submit_location,
            task_type=int(data.task_type),
            group_id=data.group_id,
            group_dict=group_dict,
        )
        params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        name = "ExpVenn" + '_'
        time_now = datetime.datetime.now()
        name += time_now.strftime("%Y%m%d_%H%M%S")
        if type(params) == dict:
            params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        main_info = dict(
            project_sn=project_sn,
            task_id=task_id,
            name=name,
            created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
            desc='Expression venn analysis main table',
            params=params,
            status="start"
        )
        main_id = self.labelfree.insert_main_table('sg_exp_venn', main_info)

        # prepare option for workflow
        if str(data.group_id).lower() == 'all':
            samples = group_dict['all']
            group_dict = OrderedDict([(x, [x]) for x in samples])
        options = {
            "exp_matrix": data.task_id,
            "group_dict": json.dumps(group_dict),
            "venn_main_id": str(main_id),
            "group": json.dumps(group_dict),
            "update_info": json.dumps({str(main_id): "sg_exp_venn"})  # to update sg_status
        }
        # prepare to file
        to_files = ["labelfree.export_exp_matrix(exp_matrix)",
                    "labelfree.export_group(group)"]

        # 把参数交给workflow运行相应的tool， 其中to_file用于准备tool的输入文件
        task_name = 'labelfree.report.exp_venn'
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
        cmd += "s/labelfree/exp_venn "
        cmd += "-b http://192.168.12.102:9090 "
        args = dict(
            task_id="labelfree",
            task_type="2",
            submit_location="expvenn",
            group_id="5aceedc15b8c915efb596de8",
            group_dict=json.dumps({"A":["A1","A2","A3"],"C":["C1","C2","C3"],"E":["E1","E2","E3"]}).replace('"', '\\"'),
            control_id="5aceee935b8c915efb596e4d",
        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)


if __name__ == '__main__':
    unittest.main()

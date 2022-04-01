# -*- coding: utf-8 -*-
import web
import json
import datetime
from mainapp.controllers.project.dia_controller import DiaController
from mbio.api.to_file.dia import *
from mainapp.libs.signature import check_sig
import unittest
import os


class ToolUpsetVennAction(DiaController):
    def __init__(self):
        super(ToolUpsetVennAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        basic_args = ["task_id", "submit_location", 'task_type']
        basic_args += ['proteinset_id', 'tool_type']
        # check arg
        for arg in basic_args:
            if not hasattr(data, arg):
                info = {'success': False, 'info': "Lack argument: {}".format(arg)}
                return json.dumps(info)
        if hasattr(data, 'submit_type'):
            submit_type = int(data.submit_type)
        else:
            submit_type = 0
        sg_task = self.dia.get_task_info(data.task_id)
        project_sn = sg_task["project_sn"]
        task_id = data.task_id
        # create main table record
        proteinset_id_list = str(data.proteinset_id).split(',')
        params = dict(
            task_id=task_id,
            submit_location=data.submit_location,
            task_type=int(data.task_type),
            proteinset_id=','.join(proteinset_id_list),
            tool_type=data.tool_type,

        )
        params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        name = "GenesetVenn" + '_'
        time_now = datetime.datetime.now()
        name += time_now.strftime("%Y%m%d_%H%M%S")
        main_info = dict(
            project_sn=project_sn,
            task_id=task_id,
            name=name,
            created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
            desc='Geneset venn analysis main table',
            params=params,
            submit_type=submit_type,
            status="start"
        )
        main_id = self.dia.insert_main_table('sg_tool_lab_upset_venn', main_info)
        # prepare option for workflow
        options = {
            "task_id": data.task_id,
            "main_id": str(main_id),
            "relate_name": name,
            "upset": data.proteinset_id,
            'params': params,
            'tool_type': data.tool_type,
            "submit_location": data.submit_location,
            "update_info": json.dumps({str(main_id): "sg_tool_lab_upset_venn"})  # to update sg_status
        }

        # 把参数交给workflow运行相应的tool， 其中to_file用于准备tool的输入文件
        to_files = ["dia_tool.export_upset_venn(upset)"]
        task_name = 'dia_v3.report.tool_upset_venn'
        self.set_sheet_data(name=task_name,
                            options=options,
                            to_file=to_files,
                            main_table_name=name,  # 设置交互分析结果目录名
                            module_type="workflow",
                            project_sn=project_sn,
                            task_id=data.task_id)

        # 运行workflow 并传回参数
        task_info = super(ToolUpsetVennAction, self).POST()
        task_info['content'] = {
            'ids': {
                'id': str(main_id),
                'name': name
                }
        }
        # 更新基因集的使用信息
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
        cmd += "s/dia_v3/tool_upset_venn "
        cmd += "-b http://wpm2.sanger.com "
        args = dict(
            task_id="t39d_v96n25npu0nq4b9nmobg5v",
            task_type="2",
            submit_location="GenesetVenn",
            proteinset_id="610b9be47a1926d23044b68c,610b9be47a1926d23044b68a,610b9be47a1926d23044b688,610b9be47a1926d23044b686",
            tool_type='upset_venn'

        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)


if __name__ == '__main__':
    unittest.main()

# -*- coding: utf-8 -*-
import web
import json
import datetime
from mainapp.controllers.project.prok_rna_controller import ProkRNAController
from mbio.api.to_file.prok_rna import *
from mainapp.libs.signature import check_sig
import unittest
import os


class ToolGoCircAction(ProkRNAController):
    def __init__(self):
        super(ToolGoCircAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        basic_args = ["task_id", "submit_location", 'task_type']
        basic_args += ['enrich_id', 'diff_id', 'cmp_group', 'tool_type']
        # check arg
        for arg in basic_args:
            if not hasattr(data, arg):
                info = {'success': False, 'info': "Lack argument: {}".format(arg)}
                return json.dumps(info)
        if hasattr(data, 'submit_type'):
            submit_type = int(data.submit_type)
        else:
            submit_type = 0
        sg_task = self.prok_rna.get_task_info(data.task_id)
        project_sn = sg_task["project_sn"]
        # create main table record
        params = dict(
            task_id=data.task_id,
            submit_location=data.submit_location,
            task_type=int(data.task_type),
            enrich_id=data.enrich_id,
            diff_id=data.diff_id,
            cmp_group=data.cmp_group,
            tool_type=data.tool_type
        )
        params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        name = "Go_Circ" + '_'
        time_now = datetime.datetime.now()
        name += time_now.strftime("%Y%m%d_%H%M%S")
        main_info = dict(
            project_sn=project_sn,
            task_id=data.task_id,
            name=name,
            created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
            desc='',
            params=params,
            submit_type=submit_type,
            status="start"
        )
        main_id = self.prok_rna.insert_main_table('sg_tool_lab_go_circ', main_info)
        # prepare option for workflow
        options = {
            "diff_exp": data.diff_id + ';' + data.cmp_group,
            "go_file": data.enrich_id,
            'tool_type': data.tool_type,
            'params': params,
            "main_id": str(main_id),
            "relate_name": name,
            'task_id': data.task_id,
            "submit_location": data.submit_location,
            "update_info": json.dumps({str(main_id): "sg_tool_lab_go_circ"})  # to update sg_status
        }

        # 把参数交给workflow运行相应的tool， 其中to_file用于准备tool的输入文件
        to_files = ["prok_rna_tool.export_go_circ_exp(diff_exp)",
                    "prok_rna_tool.export_go_circ_enrich(go_file)"]
        task_name = 'prok_rna.report.tool_go_circ'
        self.set_sheet_data(name=task_name,
                            to_file=to_files,
                            options=options,
                            main_table_name=name,  # 设置交互分析结果目录名
                            module_type="workflow",
                            project_sn=project_sn,
                            task_id=data.task_id)

        # 运行workflow 并传回参数
        task_info = super(ToolGoCircAction, self).POST()
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
        cmd += "s/prok_rna/tool_go_circ "
        cmd += "-b http://wpm2.sanger.com "
        args = dict(
            task_id="iqjq_ucq2e81kvpkr6cjpd5ifmv",
            task_type="2",
            submit_location="GoCirc",
            enrich_id="6119cf87e54aa7ba13e59323",
            diff_id="61171be634cc5373074f6620",
            cmp_group="WT|DksAppgpp",
            tool_type='Go_Circ'
        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)


if __name__ == '__main__':
    unittest.main()
# -*- coding: utf-8 -*-
import web
import json
import datetime
from mainapp.controllers.project.prok_rna_controller import ProkRNAController
from mbio.api.to_file.prok_rna import *
from mainapp.libs.signature import check_sig
import unittest
import os


class ToolGenesetPpiAction(ProkRNAController):
    def __init__(self):
        super(ToolGenesetPpiAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        basic_args = ["task_id", "submit_location", 'task_type']
        basic_args += ['geneset_id', 'tool_type']
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
        seq_path = sg_task["assemble_fa"]
        # create main table record
        params = dict(
            task_id=data.task_id,
            submit_location=data.submit_location,
            task_type=int(data.task_type),
            geneset_id=data.geneset_id,
            tool_type=data.tool_type
        )
        params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        name = "Geneset_Ppi" + '_'
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
        main_id = self.prok_rna.insert_main_table('sg_tool_lab_geneset_ppi', main_info)
        # prepare option for workflow
        options = {
            "gene_list": data.geneset_id,
            "seq_dir": os.path.dirname(seq_path),
            'tool_type': data.tool_type,
            'params': params,
            "main_id": str(main_id),
            "relate_name": name,
            'task_id': data.task_id,
            "submit_location": data.submit_location,
            "update_info": json.dumps({str(main_id): "sg_tool_lab_geneset_ppi"})  # to update sg_status
        }

        # 把参数交给workflow运行相应的tool， 其中to_file用于准备tool的输入文件
        to_files = ["prok_rna_tool.export_geneset(gene_list)"]
        task_name = 'prok_rna.report.tool_geneset_ppi'
        self.set_sheet_data(name=task_name,
                            to_file=to_files,
                            options=options,
                            main_table_name=name,  # 设置交互分析结果目录名
                            module_type="workflow",
                            project_sn=project_sn,
                            task_id=data.task_id)

        # 运行workflow 并传回参数
        task_info = super(ToolGenesetPpiAction, self).POST()
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
        cmd += "s/prok_rna/tool_geneset_ppi "
        cmd += "-b http://wpm2.sanger.com "
        args = dict(
            task_id="st1o_evgusjlkjijkgpv881df08",
            task_type="2",
            submit_location="GenesetPpi",
            geneset_id="61120f87cef2bbf53d8b4567",
            tool_type='geneset_ppi'
        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)


if __name__ == '__main__':
    unittest.main()
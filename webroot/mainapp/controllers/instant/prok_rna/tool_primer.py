# -*- coding: utf-8 -*-
import web
import json
import datetime
from mainapp.controllers.project.prok_rna_controller import ProkRNAController
from mbio.api.to_file.prok_rna import *
from mainapp.libs.signature import check_sig
import unittest
import os
import re


class ToolPrimerAction(ProkRNAController):
    def __init__(self):
        super(ToolPrimerAction, self).__init__(instant=True)

    @check_sig
    def POST(self):
        data = web.input()
        basic_args = ["task_id", "submit_location", 'task_type']
        basic_args += ['seq_id', 'updownstream', 'tool_type']
        # check arg
        for arg in basic_args:
            if not hasattr(data, arg):
                info = {'success': False, 'info': "Lack argument: {}".format(arg)}
                return json.dumps(info)
        if hasattr(data, 'submit_type'):
            submit_type = int(data.submit_type)
        else:
            submit_type = 0
        # sg_task = self.prok_rna.get_task_info(data.task_id)
        # project_sn = sg_task["project_sn"]
        # task_id = data.task_id
        sg_task = self.prok_rna.get_task_info(data.task_id)
        project_sn = sg_task["project_sn"]
        ref_genome = sg_task['ref_genome']

        rock_info = sg_task['rock_index']

        if re.match(r'^\w+://\S+/.+$', rock_info):
            inter_dir = self.create_tmp_dir(data.task_id, "rock_info/")
            rock_info_fna = self.download_from_s3(os.path.join(rock_info, 'genome_fna.db.sqlite3'), inter_dir=inter_dir)
            rock_info_bed = self.download_from_s3(os.path.join(rock_info, 'ptt.bed'), inter_dir=inter_dir)
            rock_info = rock_info_bed.split('ptt.bed')[0]

        ptt_path = os.path.join(rock_info, 'ptt.bed')
        genome_path = os.path.join(rock_info, 'genome_fna.db.sqlite3')
        # create main table record
        params = dict(
            task_id=data.task_id,
            submit_location=data.submit_location,
            task_type=int(data.task_type),
            seq_id = data.seq_id,
            updownstream = data.updownstream,
            tool_type = data.tool_type
        )
        params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        name = "Seqdown" + '_'
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
        main_id = self.prok_rna.insert_main_table('sg_tool_lab_seqdown', main_info)
        # prepare option for workflow
        options = {
            "task_id": data.task_id,
            'params': params,
            'tool_type': data.tool_type,
            # "submit_location": data.submit_location,
            # 'update_info': json.dumps(update_info),
            "main_id": str(main_id),
            "seq_id": data.seq_id,
            "upstream": int(data.updownstream),
            "downstream": int(data.updownstream),
            "ref_genome": ref_genome,
            "ptt_path": ptt_path,
            "genome_path": genome_path,
            "submit_location": data.submit_location,
            "update_info": json.dumps({str(main_id): "sg_tool_lab_seqdown"})  # to update sg_status
        }

        # 把参数交给workflow运行相应的tool， 其中to_file用于准备tool的输入文件
        task_name = 'prok_rna.report.tool_primer'
        self.set_sheet_data(name=task_name,
                            options=options,
                            main_table_name=name,  # 设置交互分析结果目录名
                            module_type="workflow",
                            project_sn=project_sn,
                            task_id=data.task_id)

        # 运行workflow 并传回参数
        task_info = super(ToolPrimerAction, self).POST()
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
        cmd = 'python /mnt/ilustre/users/sanger-dev/biocluster/bin/webapitest.py '
        cmd += 'post '
        cmd += "-fr no "
        cmd += '-c {} '.format("client03")
        cmd += '-dbversion {} '.format(0)
        cmd += "s/prok_rna/tool_primer "
        cmd += "-b http://bcl.tsg.com "
        args = dict(
            task_id="tsg_248821",
            task_type="2",
            submit_location="seqdown",
            seq_id="YEP0002",
            tool_type='primer',
            updownstream = '100'

        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)


if __name__ == '__main__':
    unittest.main()

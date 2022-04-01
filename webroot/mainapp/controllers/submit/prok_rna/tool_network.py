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


class ToolNetworkAction(ProkRNAController):
    def __init__(self):
        super(ToolNetworkAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        basic_args = ["task_id", "submit_location", 'task_type']
        basic_args += ['group_dict', 'diff_id', 'exp_id', 'cmp_list', 'gene_type', 'significant', 'regulate', 'stat_type', 'stat_type_operation',
                       'stat_type_value', 'log2fc', 'log2fc_operation', 'log2fc_value', 'geneset_id', 'tool_type']
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
            group_dict=data.group_dict,
            diff_id=data.diff_id,
            exp_id=data.exp_id,
            cmp_list=data.cmp_list,
            gene_type=data.gene_type,
            significant=data.significant,
            regulate=data.regulate,
            stat_type=data.stat_type,
            stat_type_operation=data.stat_type_operation,
            stat_type_value=data.stat_type_value,
            log2fc=data.log2fc,
            log2fc_operation=data.log2fc_operation,
            log2fc_value=data.log2fc_value,
            geneset_id=data.geneset_id,
            tool_type=data.tool_type
        )
        params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        name = "Network" + '_'
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
        main_id = self.prok_rna.insert_main_table('sg_tool_lab_network', main_info)
        # prepare option for workflow
        options = {
            "group_dict": data.group_dict,
            "network_file": data.diff_id,
            "exp_id": data.exp_id,
            'cmp_list': data.cmp_list,
            'gene_type': data.gene_type,
            'significant': data.significant,
            'regulate': data.regulate,
            'stat_type': data.stat_type,
            'stat_type_operation': data.stat_type_operation,
            'stat_type_value': data.stat_type_value,
            'log2fc': data.log2fc,
            'log2fc_operation': data.log2fc_operation,
            'log2fc_value': data.log2fc_value,
            'geneset_id': data.geneset_id,
            'tool_type': data.tool_type,
            'params': params,
            "main_id": str(main_id),
            "relate_name": name,
            'task_id': data.task_id,
            "submit_location": data.submit_location,
            "update_info": json.dumps({str(main_id): "sg_tool_lab_network"})  # to update sg_status
        }

        # 把参数交给workflow运行相应的tool， 其中to_file用于准备tool的输入文件
        to_files = ["prok_rna_tool.export_network_diff(network_file)"]
        task_name = 'prok_rna.report.tool_network'
        self.set_sheet_data(name=task_name,
                            to_file=to_files,
                            options=options,
                            main_table_name=name,  # 设置交互分析结果目录名
                            module_type="workflow",
                            project_sn=project_sn,
                            task_id=data.task_id)

        # 运行workflow 并传回参数
        task_info = super(ToolNetworkAction, self).POST()
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
        cmd += '-dbversion {} '.format(1)
        cmd += "s/prok_rna/tool_network "
        cmd += "-b http://bcl.tsg.com "
        args = dict(
            task_id="gamk_3vf5dqd4o2og8rd9cm7pq0",
            task_type="2",
            submit_location="network",
            group_dict=r'{"E1": ["E1_1", "E1_2", "E1_3"], "N1": ["N1_1", "N1_2", "N1_3"]}'.replace('"', '\\"'),
            diff_id='6058e89717b2bf47ca3029da',
            exp_id="6058e88317b2bf47ca2ff9be",
            cmp_list='E1|N1',
            gene_type='mRNA',
            significant='',
            regulate='',
            stat_type='pvalue',
            stat_type_operation='',
            stat_type_value='',
            log2fc='fc',
            log2fc_operation='',
            log2fc_value='',
            geneset_id='',
            tool_type='network'

        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)


if __name__ == '__main__':
    unittest.main()

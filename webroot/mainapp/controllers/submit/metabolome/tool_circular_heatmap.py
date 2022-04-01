# -*- coding: utf-8 -*-
import web
import json
import datetime
from mainapp.controllers.project.metabolome_controller import MetabolomeController
from mbio.api.to_file.metabolome import *
from mainapp.libs.signature import check_sig
from mainapp.models.mongo.metabolome import Metabolome
import unittest
import os
from collections import OrderedDict


class ToolCircularHeatmapAction(MetabolomeController):
    def __init__(self):
        super(ToolCircularHeatmapAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        basic_args = ["task_id", "submit_location", 'task_type']
        basic_args += ['exp_id', 'group_dict', 'metabset_id', 'tool_type']
        # check arg
        for arg in basic_args:
            if not hasattr(data, arg):
                info = {'success': False, 'info': "Lack argument: {}".format(arg)}
                return json.dumps(info)
        if hasattr(data, 'submit_type'):
            submit_type = int(data.submit_type)
        else:
            submit_type = 0
        metabolome = Metabolome()
        project_sn, ms_type, exp_name = metabolome.get_metab_info(data.exp_id, data.task_id)

        if hasattr(data, "table_type") and exp_name == "raw":
            if data.table_type == "mix":
                info = {"success": False, "info": "原始表没有合并表,请选择正确参数!", 'code':'C2300801'}
                return json.dumps(info)

        # create main table record
        params = dict(
            task_id=data.task_id,
            submit_location=data.submit_location,
            task_type=int(data.task_type),
            group_dict=data.group_dict,
            exp_id=data.exp_id,
            metabset_id=data.metabset_id,
            tool_type=data.tool_type
        )
        if ms_type == "LC":
            if hasattr(data, "table_type"):
                # params["table_type"] = data.table_type
                table_type = data.table_type
            else:
                table_type = 'mix'
        elif ms_type == "GC":
            # params["table_type"] = 'pos'
            table_type = 'pos'

        name = "Circular_Heatmap" + '_'
        time_now = datetime.datetime.now()
        name += time_now.strftime("%Y%m%d_%H%M%S")
        if hasattr(data, 'table_name'):
            table_name = data.table_name
            params.update({'table_name': data.table_name})
        else:
            table_name = name
        params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        main_info = dict(
            project_sn=project_sn,
            task_id=data.task_id,
            name=table_name,
            created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
            desc='',
            params=params,
            submit_type=submit_type,
            status="start"
        )
        main_id = metabolome.insert_main_table('sg_tool_lab_circular_heatmap', main_info)
        # prepare option for workflow
        options = {
            "group_dict": data.group_dict,
            "metab_table": data.exp_id + ';' + data.metabset_id + ';' + table_type,
            'tool_type': data.tool_type,
            'params': params,
            "main_id": str(main_id),
            "relate_name": table_name,
            'task_id': data.task_id,
            "submit_location": data.submit_location,
            "ms_type": ms_type,
            "update_info": json.dumps({str(main_id): "sg_tool_lab_circular_heatmap"})  # to update sg_status
        }

        # 把参数交给workflow运行相应的tool， 其中to_file用于准备tool的输入文件
        to_files = ["metabolome_tool.export_metabset_exp(metab_table)"]
        task_name = 'metabolome.report.tool_circular_heatmap'
        self.set_sheet_data(name=task_name,
                            to_file=to_files,
                            options=options,
                            main_table_name=table_name,  # 设置交互分析结果目录名
                            module_type="workflow",
                            project_sn=project_sn,
                            task_id=data.task_id)

        # 运行workflow 并传回参数
        task_info = super(ToolCircularHeatmapAction, self).POST()
        task_info['content'] = {
            'ids': {
                'id': str(main_id),
                'name': table_name
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
        cmd += "s/metabolome/tool_circular_heatmap "
        cmd += "-b http://wpm2.sanger.com "
        args = dict(
            task_id="kfns_pm0lnh9v8267nae9crjduv",
            task_type="2",
            submit_location="CircularHeatmap",
            group_dict=r'{"A": ["A1", "A2", "A3"], "R": ["R1", "R2", "R3"]}'.replace('"', '\\"'),
            metabset_id='611377ff5a64f865538b9ef3',
            tool_type='circular_heatmap',
            exp_id='611377ff5a64f865538b9e5b',
        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)


if __name__ == '__main__':
    unittest.main()
#!/usr/bin/python
# -*- coding: utf-8 -*-
# __author__ = 'konghualei,20170426'

import web
import json
from mainapp.libs.signature import check_sig
import datetime
from mainapp.controllers.project.tool_lab_controller import ToolLabController
from mainapp.models.mongo.tool_lab import *
import os
import unittest
from collections import OrderedDict

class TableKitStandardAction(ToolLabController):
    def __init__(self):
        super(TableKitStandardAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        params_name = ['task_id', 'table', 'sep', 'observe', 'feature', 'submit_location']
        for param in params_name:
            if not hasattr(data, param):
                info = {'success': False, 'info': "缺少%s参数！" % param, 'variables': params_name}
                return json.dumps(info)
        # inter_dir = self.create_tmp_dir(data.task_id, 'table_kit_standard')
        # table = self.download_from_s3(data.table, inter_dir=inter_dir)

        params_json = {
            'table': data.table,
            # 'manipulating': data.manipulating,
            'sep': data.sep,
            'observe': data.observe,
            'feature': data.feature,
            'task_id':data.task_id,
            'submit_location' : data.submit_location
        }

        if hasattr(data, "guard"):
            params_json.update({
                "guard": data.guard
            })

        params = json.dumps(params_json, sort_keys=True, separators=(',', ':'))

        name = "Standard" + '_'
        time_now = datetime.datetime.now()
        name += time_now.strftime("%Y%m%d_%H%M%S")
        main_info = dict(
            project_sn='tool',
            task_id=data.task_id,
            name=name,
            # version="v1",
            created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
            desc='Table standard main table',
            params=params,
            status="start",

        )
        print main_info
        main_id = self.tool_lab.insert_main_table('Table_standard', main_info)
        options = {
            'table': data.table,
            # 'manipulating': data.manipulating,
            'sep': data.sep,
            'observe': data.observe,
            'feature': data.feature,
            # 'task_id': data.task_id,
            'update_info' : json.dumps({str(main_id): "table_kit_standard"}),#sg_status
            'standard_id' : str(main_id)

        }

        task_name = 'tool_lab.table_kit_standard'
        self.set_sheet_data(name=task_name,
                            options=options,
                            main_table_name=name,  # 设置交互分析结果目录名
                            module_type="workflow",
                            to_file=None,
                            project_sn='tool',
                            task_id=data.task_id)
        task_info = super(TableKitStandardAction, self).POST()
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
        cmd += "s/tool_lab/table_kit_standard "
        cmd += "-b http://bcl.tsg.com "
        args = dict(
            task_id="Standard",
            task_type="2",
            submit_location="Standard_upload",
            # manipulating='standard',
            sep='tab',
            observe='sum',
            feature='standard_scale',
            table='/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/table.txt',

        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)


if __name__ == '__main__':
    unittest.main()
# -*- coding: utf-8 -*-
# __author__ = fiona
# last modified by shicaiping at 20180515

import re, os, Bio, argparse, sys, fileinput
from mainapp.libs.signature import check_sig
from mainapp.libs.param_pack import *
from mbio.api.to_file.ref_rna_v2 import *
from mainapp.controllers.project.ref_rna_v2_controller import RefRnaV2Controller
import web
import subprocess
import json
import datetime
import unittest


class RmatsModelAction(RefRnaV2Controller):
    def __init__(self):
        super(RmatsModelAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        basic_args = ["task_id", "submit_location", 'task_type']
        basic_args += ["splicing_id", "gene_id", "as_type", "event_type"]
        for arg in basic_args:
            if not hasattr(data, arg):
                variables = []
                variables.append(arg)
                info = {'success': False, 'info': '%s参数缺少!' % arg, 'code': 'C2902401', 'variables': variables}
                return json.dumps(info)
        event_type = data.event_type
        if event_type not in ["A3SS", "A5SS", "MXE", "SE", "RI"]:
            variables = []
            variables.append(event_type)
            info = {'success': False, 'info': "%s分析不存在" % event_type, 'code': 'C2902402', 'variables': variables}
            return json.dumps(info)
        splicing_rmats_info = self.ref_rna_v2.get_main_info(data.splicing_id, 'sg_splicing_rmats', data.task_id)
        if not splicing_rmats_info:
            info = {"success": False, "info": "splicing_rmats主表信息不存在，请检查参数是否正确！", 'code': 'C2902403', 'variables': ''}
            return json.dumps(info)
        try:
            rmats_out_root_dir = splicing_rmats_info['result_dir']
        except Exception as e:
            info = {"success": False, "info": "splicing_rmats主表信息中不存在结果目录信息！", 'code': 'C2902404', 'variables': ''}
            return json.dumps(info)
        if not rmats_out_root_dir.endswith("/"):
            rmats_out_root_dir += "/"
        #task_info = self.ref_rna_v2.get_task_info(splicing_rmats_info['task_id'])
        task_type = int(data.task_type)
        task_name = 'ref_rna_v2.report.rmats_model'
        my_param = {'splicing_id': data.splicing_id, 'gene_id': data.gene_id,
                    'event_type': data.event_type, 'as_type': data.as_type, 'task_id': data.task_id,
                    'task_type': task_type, 'submit_location': data.submit_location}
        main_table_name =  "SplicingRmatsModel_" + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S"))
        task_id = splicing_rmats_info['task_id']
        task_info = self.ref_rna_v2.get_task_info(data.task_id)
        project_sn = task_info['project_sn']
        mongo_data = dict([
            ('task_id', data.task_id),
            ('status', 'start'),
            ('name', main_table_name),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ("params", json.dumps(my_param, sort_keys=True, separators=(',', ':'))),
            ("graph_dir", ""),
            ("png_files", ""),
            ("pdf_files", ""),
        ])
        collection_name = "sg_splicing_rmats_model"
        self.ref_rna_v2.delete_main_table(collection_name, data.task_id)
        main_table_id = self.ref_rna_v2.insert_main_table(collection_name, mongo_data)
        new_task_id = self.ref_rna_v2.get_new_id(data.task_id)
        main_table_data = {'run_id': new_task_id}
        update_info = {str(main_table_id): collection_name}
        update_info = json.dumps(update_info)
        options = {
            "splicing_id": data.splicing_id,
            "main_table_data": main_table_data,
            "update_info": update_info,
            "gene_id": data.gene_id,
            "event_type": data.event_type,
            "as_type": data.as_type,
            "rmats_out_root_dir": self.use_s3(rmats_out_root_dir),
            "task_id": data.task_id,
        }
        self.set_sheet_data(name=task_name, options=options, main_table_name=main_table_name, task_id=task_id,
                            project_sn=project_sn, module_type='workflow', new_task_id=new_task_id)
        task_info = super(RmatsModelAction, self).POST()
        task_info['content'] = {'ids': {'id': str(main_table_id), 'name': main_table_name}}
        if 'group_id' in data and str(data.group_id).lower() != 'all':
            _ = self.ref_rna_v2.update_group_is_use(data.task_id, data.group_id)
        if 'control_id' in data:
            _ = self.ref_rna_v2.update_group_compare_is_use(data.task_id, data.control_id)
        return json.dumps(task_info)

class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test_this(self):
        import os
        cmd = 'python /mnt/ilustre/users/sanger-dev/biocluster/bin/webapitest.py '
        cmd += 'post '
        cmd += "-fr no "
        cmd += '-c {} '.format("client03")
        cmd += "s/ref_rna_v2/rmats_model "
        cmd += "-b http://192.168.12.102:9090 "
        args = dict(
            task_id="RefrnaV2_7320",
            submit_location="rmatsmodel",
            task_type="2",
            splicing_id="5b1629dea4e1af36ff8be215",
            gene_id="ENSMUSG00000018900",
            as_type="JC",
            event_type="SE",
        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)

if __name__ == '__main__':
    unittest.main()

# -*- coding: utf-8 -*-
# __author__ = shicaiping

import web
import json
from mainapp.libs.signature import check_sig
from bson.objectid import ObjectId
from mainapp.libs.param_pack import *
from mbio.api.to_file.ref_rna_v2 import *
from mainapp.models.mongo.ref_rna_v2 import RefRnaV2
from mainapp.controllers.project.ref_rna_v2_controller import RefRnaV2Controller
from biocluster.config import Config
from mainapp.libs.param_pack import group_detail_sort
from collections import OrderedDict
import datetime
import unittest
import os,re

class SplicingRmatsAction(RefRnaV2Controller):
    '''
    可变剪切接口
    '''
    def __init__(self):
        super(SplicingRmatsAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        basic_args = ["task_id", "submit_location", 'task_type']
        basic_args += ['group_id', 'group_dict', "compare_plan","control_id"]
        for arg in basic_args:
            if not hasattr(data, arg):
                variables = []
                variables.append(arg)
                info = {'success': False, 'info': '%s参数缺少!' % arg, 'code': 'C2902701', 'variables': variables}
                return json.dumps(info)
        task_name = 'ref_rna_v2.report.rmats'
        case_group_name = data.compare_plan.split('|')[0]
        print(case_group_name)
        control_group_name = data.compare_plan.split('|')[1]
        print(control_group_name)
        group_dict = json.loads(data.group_dict, object_pairs_hook=OrderedDict)
        print(group_dict)
        case_group_members = group_dict[case_group_name]
        print(case_group_members)
        control_group_members = group_dict[control_group_name]
        print(control_group_members)

        case_group_bam_lst = []
        control_group_bam_lst = []
        for case_sp in sorted(case_group_members):
            bam_path = RefRnaV2().get_bam_path(case_sp, data.task_id)
            case_group_bam_lst.append(bam_path)
        for control_sp in sorted(control_group_members):
            bam_path = RefRnaV2().get_bam_path(control_sp, data.task_id)
            control_group_bam_lst.append(bam_path)
        print(case_group_bam_lst)
        print(control_group_bam_lst)
        case_group_bam_str = ','.join([p.strip() for p in case_group_bam_lst])
        control_group_bam_str = ','.join([p.strip() for p in control_group_bam_lst])
        my_param = {}
        my_param['task_id'] = data.task_id
        my_param['group_dict'] = group_detail_sort(data.group_dict)
        my_param['group_id'] = data.group_id
        my_param['compare_plan'] = data.compare_plan
        my_param['control_id'] = data.control_id
        my_param['submit_location'] = data.submit_location
        my_param['task_type'] = int(data.task_type)
        task_info = self.ref_rna_v2.get_task_info(data.task_id)
        if task_info:
            main_table_name = "SplicingRmats_" + case_group_name + "_vs_" + control_group_name + "_" + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S"))
            task_id = data.task_id
            project_sn = task_info["project_sn"]
            ref_gtf = task_info["ref_gtf"]
            app_dir = Config().SOFTWARE_DIR
            if not os.path.exists(ref_gtf):
                ref_gtf = os.path.join(app_dir, ref_gtf.split("/app/")[1])
            group = {case_group_name: case_group_members, control_group_name: control_group_members}
            mongo_data = [
                ('task_id', task_info['task_id']),
                ('status', 'start'),
                ('desc', "可变剪接rmats计算主表"),
                ('name', main_table_name),
                ('compare_plan', data.compare_plan),
                ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
                ("params", json.dumps(my_param, sort_keys=True, separators=(',', ':')))
            ]
            collection_name = "sg_splicing_rmats"
            main_table_id = self.ref_rna_v2.insert_main_table(collection_name, mongo_data)
            update_info = {str(main_table_id): collection_name}
            update_info = json.dumps(update_info)
            new_task_id = self.ref_rna_v2.get_new_id(data.task_id)
            main_table_data = {'run_id': new_task_id}
            options = {
                'main_table_data': main_table_data,
                "update_info": update_info,
                "ref_gtf": self.use_s3(ref_gtf),
                'control_file': data.control_id,
                "group_id": data.group_id,
                "group_dict": data.group_dict,
                "splicing_id": str(main_table_id)
            }
            if data.group_id != 'all':
                options.update({
                    "case_group_bam_str": case_group_bam_str,
                    "control_group_bam_str": control_group_bam_str,
                    "case_group_name": case_group_name,
                    "control_group_name": control_group_name
                })
            else:
                info = {"success": False, "info": "不能选单个样本作为一个组 组必须包含多个样本", 'code': 'C2902702', 'variables': ''}
                return json.dumps(info)
            self.set_sheet_data(name=task_name, options=options, main_table_name=main_table_name, task_id=task_id,
                                project_sn=project_sn, module_type='workflow', params=my_param, new_task_id=new_task_id)
            task_info = super(SplicingRmatsAction, self).POST()
            task_info['content'] = {'ids': {'id': str(main_table_id), 'name': main_table_name}}
            if 'group_id' in data and str(data.group_id).lower() != 'all':
                _ = self.ref_rna_v2.update_group_is_use(data.task_id, data.group_id)
            if 'control_id' in data:
                _ = self.ref_rna_v2.update_group_compare_is_use(data.task_id, data.control_id)
            return json.dumps(task_info)
        else:
            info = {"success": False, "info": "sg_task表不存在，请确认参数是否正确！!", 'code': 'C2902703', 'variables': ''}
            return json.dumps(info)

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
        cmd += "s/ref_rna_v2/splicing_rmats "
        cmd += "-b http://192.168.12.102:9090 "
        args = dict(
            task_id="RefrnaV2_7320",
            submit_location="splicingrmats",
            task_type="2",
            group_id="5afd30bba4e1af301da5f379",
            control_id="5afd30bba4e1af301da5f37a",
            group_dict=json.dumps({"A1":["A1_1", "A1_2", "A1_3"],"A2":["A2_1", "A2_2", "A2_3"]}).replace('"', '\\"'),
            compare_plan="A1|A2",
        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)

if __name__ == '__main__':
    unittest.main()

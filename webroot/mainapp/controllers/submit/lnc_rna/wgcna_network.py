# -*- coding: utf-8 -*-
import web
import json
import datetime
from collections import OrderedDict
from mainapp.controllers.project.lnc_rna_controller import LncRnaController
from mbio.api.to_file.lnc_rna import *
from mainapp.libs.signature import check_sig
import unittest
import os
import types
from bson.objectid import ObjectId

class WgcnaNetworkAction(LncRnaController):
    def __init__(self):
        super(WgcnaNetworkAction, self).__init__(instant=False)
        self.task_name = 'lnc_rna.report.wgcna_network'
        self.expected_args = ["task_id", "submit_location", 'task_type']
        self.expected_args += ['module', 'wgcna_relate_id', 'threshold', 'top', 'exp_level']
        self.input_data = web.input()

    @check_sig
    def POST(self):
        check = self.check_params()
        if check is not True:
            return check
        packed_params = self.pack_params()
        main_id, main_table_name = self.create_main_table(packed_params, "sg_wgcna_network")
        self.prepare_workflow_options(main_id, "sg_wgcna_network", main_table_name)
        task_info = super(WgcnaNetworkAction,self).POST()
        task_info['content'] = {'ids': {'id': str(main_id), 'name': main_table_name}}
        if 'group_id' in self.input_data and str(self.input_data.group_id).lower() != 'all':
            _ = self.lnc_rna.update_group_is_use(self.input_data.task_id, self.input_data.group_id)
        if 'control_id' in self.input_data:
            _ = self.lnc_rna.update_group_compare_is_use(self.input_data.task_id, self.input_data.control_id)
        return json.dumps(task_info)

    def check_params(self):
        for arg in self.expected_args:
            if not hasattr(self.input_data, arg):
                variables = []
                variables.append(arg)
                info = {'success': False, 'info': "Lack argument: %s" % arg, 'code': 'C2903301', 'variables': variables}
                return json.dumps(info)
        # do some other check
        if len(self.input_data.module) < 1:
            info = {'success': False, 'info': "module: 选择的模块为空", 'code': 'C2903302', 'variables': ''}
            return json.dumps(info)
        return True

    def pack_params(self):
        params = dict(
            task_id=self.input_data.task_id,
            submit_location=self.input_data.submit_location,
            task_type=int(self.input_data.task_type),
            wgcna_relate_id=self.input_data.wgcna_relate_id,
            threshold=self.input_data.threshold,
            module=self.input_data.module,
            top=self.input_data.top,
            exp_level=self.input_data.exp_level,
        )
        return json.dumps(params, sort_keys=True, separators=(',',':'))

    def create_main_table(self, packed_params, db_name):
        result_info = self.lnc_rna.get_main_info(self.input_data.wgcna_relate_id, 'sg_wgcna_relate', self.input_data.task_id)
        self.project_sn = result_info["project_sn"]
        self.wgcna_module_id = result_info['wgcna_module_id']
        # name = "WgcnaNetwork" + '_' + self.input_data.module.replace(",", "_") + '_'
        # 名称太长导致页面对象存储下载报错
        name = "WgcnaNetwork" + '_'
        time_now = datetime.datetime.now()
        name += time_now.strftime("%Y%m%d_%H%M%S")
        if isinstance(self.input_data.wgcna_relate_id, types.StringTypes):
            wgcna_relate_id = ObjectId(self.input_data.wgcna_relate_id)
        elif isinstance(self.input_data.wgcna_relate_id, ObjectId):
            wgcna_relate_id = self.input_data.wgcna_relate_id
        else:
            raise Exception("wgcna_prepare_id参数必须为字符串或者ObjectId类型!")
        main_info = dict(
            project_sn=self.project_sn,
            task_id=self.input_data.task_id,
            name=name,
            created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
            wgcna_module_id=self.wgcna_module_id,
            wgcna_relate_id=wgcna_relate_id,
            desc='wgcna network analysis',
            params=packed_params,
            status="start")
        main_id = self.lnc_rna.insert_main_table(db_name, main_info)
        return main_id, name

    def prepare_workflow_options(self, main_id, db_name, main_table_name):
        relate_output = self.lnc_rna.get_main_info(self.input_data.wgcna_relate_id,'sg_wgcna_relate', self.input_data.task_id)['output_dir']
        if not relate_output.endswith('/'):
            relate_output += '/'
        module_output = self.lnc_rna.get_main_info(self.wgcna_module_id, 'sg_wgcna_module', self.input_data.task_id)['output_dir']
        if not module_output.endswith('/'):
            module_output += '/'
        options = {
            "main_id": str(main_id),
            "module": self.input_data.module,
            "threshold": self.input_data.threshold,
            "top": self.input_data.top,
            "step3output":relate_output,
            "step2output":module_output,
            "update_info": json.dumps({str(main_id): db_name})  # to update sg_status
        }
        # prepare to file
        to_files = []
        # 把参数交给workflow运行相应的tool， 其中to_file用于准备tool的输入文件
        self.set_sheet_data(
            name=self.task_name,
            options=options,
            main_table_name=main_table_name,  # 设置交互分析结果目录名
            module_type="workflow",
            to_file=to_files,
            project_sn=self.project_sn,
            task_id=self.input_data.task_id
        )


class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test_this(self):
        cmd = 'python /mnt/ilustre/users/sanger-dev/biocluster/bin/webapitest.py '
        cmd += 'post '
        cmd += "-fr no "
        cmd += '-c {} '.format("client03")
        cmd += "s/lnc_rna/wgcna_network "
        cmd += "-b http://192.168.12.101:9090 "
        args = dict(
            task_id="lnc_rna",
            task_type="2",
            submit_location="WgcnaNetwork",
            wgcna_relate_id="5ca163ff17b2bf1b708b2d37",
            module='grey,blue',
            threshold='0.03',
            exp_level="gene",
            top='30',
        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)


if __name__ == '__main__':
    unittest.main()

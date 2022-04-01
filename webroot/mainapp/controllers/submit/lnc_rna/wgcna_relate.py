# -*- coding: utf-8 -*-
import web
import json
import datetime
from collections import OrderedDict
from mainapp.controllers.project.lnc_rna_controller import LncRnaController
from mbio.api.to_file.lnc_rna import *
from mainapp.libs.signature import check_sig
import unittest
import os,re
import chardet
import types
from bson.objectid import ObjectId

class WgcnaRelateAction(LncRnaController):
    def __init__(self):
        super(WgcnaRelateAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        basic_args = ["task_id", "submit_location", 'task_type']
        basic_args += ['exp_level', 'wgcna_module_id', 'trait_path', 'corr_method', 'trait_type', "file_id"]
        # check arg
        for arg in basic_args:
            if not hasattr(data, arg):
                variables = []
                variables.append(arg)
                info = {'success': False, 'info': "Lack argument: %s" % arg, 'code': 'C2903601', 'variables': variables}
                return json.dumps(info)
        main_info = self.lnc_rna.get_main_info(data.wgcna_module_id, 'sg_wgcna_module', data.task_id)
        project_sn = main_info["project_sn"]
        task_id = data.task_id
        # get abs path of input file and check
        client = data.client if hasattr(data,"client") else web.ctx.env.get('HTTP_CLIENT')
        if client == 'client01':
            target_dir = 'data'
        else:
            target_dir = 'tsanger-data'
        base_path = "/mnt/ilustre/{}/".format(target_dir)
        trait_path = base_path + data.trait_path
        if os.path.exists(trait_path):
            pass
        elif re.match(r'^\w+://\S+/.+$', data.trait_path):
            inter_dir = self.create_tmp_dir(data.task_id, "wgcna_pipeline/")
            trait_path = self.download_from_s3(data.trait_path, inter_dir=inter_dir)
        else:
            raise "文件传递格式错误 {}".format(data.trait_path)

        with open(trait_path, 'rb') as f:
            def is_number(s):
                try:
                    float(s)
                    return True
                except ValueError:
                    pass
            # check content
            header = f.readline()
            encoding_format = chardet.detect(header)
            if not (encoding_format['confidence'] == 1 and
                    encoding_format['encoding'].lower().startswith(('ascii', 'utf-8'))):
                return json.dumps({'success': False, 'info': "文件编码格式必须是ascii或utf-8", 'code': 'C2903602', 'variables': ''})
            first_sample, test_one = f.readline().strip().split()[0:2]
            if data.trait_type.lower() == "discrete":
                if is_number(test_one):
                    return json.dumps({'success': False, 'info': "选择非连续类型时，表型数据不应该是数字", 'code': 'C2903603', 'variables': ''})
            else:
                if not is_number(test_one):
                    print(trait_path)
                    return json.dumps({'success': False,'info': "选择连续类型时，表型数据必须是数字", 'code': 'C2903604', 'variables': ''})
            # check sample number
            trait_samples = [first_sample] + [x.strip().split()[0] for x in f if x.strip()]
            if set(main_info['samples']) != set(trait_samples):
                return json.dumps({'success': False, 'info': "表型样本与预处理样本不一致！", 'code': 'C2903605', 'variables': ''})

        # create main table record
        params = dict(
            task_id=task_id,
            submit_location=data.submit_location,
            task_type=int(data.task_type),
            wgcna_module_id=data.wgcna_module_id,
            exp_level=data.exp_level,
            corr_method=data.corr_method,
            trait_type=data.trait_type,
            trait_path=data.trait_path,
            file_id=data.file_id,
        )
        params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        name = "WgcnaRelate" + '_' + data.exp_level[0].upper() + '_'
        time_now = datetime.datetime.now()
        name += time_now.strftime("%Y%m%d_%H%M%S")
        if isinstance(data.wgcna_module_id, types.StringTypes):
            wgcna_module_id = ObjectId(data.wgcna_module_id)
        elif isinstance(data.wgcna_module_id, ObjectId):
            wgcna_module_id = data.wgcna_module_id
        else:
            raise Exception("wgcna_prepare_id参数必须为字符串或者ObjectId类型!")
        module_info = self.lnc_rna.get_main_info_by_record('sg_wgcna_module', main_id=wgcna_module_id, task_id=data.task_id)
        module_rdata = os.path.join(module_info["output_dir"], 'blockwiseModules_result.RData')
        main_info = dict(
            project_sn=project_sn,
            task_id=task_id,
            name=name,
            created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
            wgcna_module_id=wgcna_module_id,
            desc='wgcna relation analysis main table',
            type=data.exp_level,
            params=params,
            status="start"
        )
        main_id = self.lnc_rna.insert_main_table('sg_wgcna_relate', main_info)

        # prepare option for workflow
        options = {
            "exp_eigengenes": data.wgcna_module_id,
            "main_id": str(main_id),
            "trait_path": trait_path,
            "trait_type": data.trait_type,
            "corr_method": data.corr_method,
            "block_Rdata": module_rdata,
            "exp_level": data.exp_level,
            "update_info": json.dumps({str(main_id): "sg_wgcna_relate"})  # to update sg_status
        }
        # prepare to file
        to_files = ["lnc_rna.export_wgcna_relate_input(exp_eigengenes)"]

        # 把参数交给workflow运行相应的tool， 其中to_file用于准备tool的输入文件
        task_name = 'lnc_rna.report.wgcna_relate'
        self.set_sheet_data(name=task_name,
                            options=options,
                            main_table_name=name,  # 设置交互分析结果目录名
                            module_type="workflow",
                            to_file=to_files,
                            project_sn=project_sn,
                            task_id=data.task_id)

        # 运行workflow 并传回参数
        task_info = super(WgcnaRelateAction, self).POST()
        task_info['content'] = {
            'ids': {
                'id': str(main_id),
                'name': name
                }
        }
        # task_info['group_dict'] = group_dict
        if 'group_id' in data and str(data.group_id).lower() != 'all':
            _ = self.lnc_rna.update_group_is_use(data.task_id, data.group_id)
        if 'control_id' in data:
            _ = self.lnc_rna.update_group_compare_is_use(data.task_id, data.control_id)
        return json.dumps(task_info)


class TestFunction(unittest.TestCase):
    """
    This is test for the workflow. Just run this script to do test.
    """
    def test_this(self):
        cmd = 'python /mnt/ilustre/users/sanger-dev/biocluster/bin/webapitest.py '
        cmd += 'post '
        cmd += "-fr no "
        cmd += '-c {} '.format("client03")
        cmd += "s/lnc_rna/wgcna_relate "
        cmd += "-b http://192.168.12.102:9090 "
        args = dict(
            task_id="tsg_28226",
            task_type="2",
            submit_location="wgcna_relate",
            wgcna_module_id="5abb54b2a4e1af3d91b8eb3e",
            exp_level="G",
            trait_path="/mnt/ilustre/users/sanger-dev/sg-users/deqing/projects_test_input/traits.xls",
            trait_type="discrete",
            corr_method="pearson",
        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)


if __name__ == '__main__':
    unittest.main()

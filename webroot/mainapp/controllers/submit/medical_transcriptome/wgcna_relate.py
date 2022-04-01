# -*- coding: utf-8 -*-

import datetime
import json
import os
import re
import types
import unittest
import web
from bson.objectid import ObjectId
from mainapp.controllers.project.medical_transcriptome_controller import MedicalTranscriptomeController
from mainapp.libs.signature import check_sig
from mbio.packages.rna.table_check import Table
from mbio.api.to_file.medical_transcriptome import *


class WgcnaRelateAction(MedicalTranscriptomeController):
    def __init__(self):
        super(WgcnaRelateAction, self).__init__(instant=False)

    def clean_trait(self, trait_path):
        with open(trait_path, 'rb') as f, open(trait_path + ".checked.txt", 'w') as fo:
            header = f.readline()
            cols = header.strip().lstrip("#").split()
            fo.write("\t".join(cols) + "\n")
            for line in f:
                cols = line.strip().split()
                fo.write("\t".join(cols) + "\n")
        return trait_path + ".checked.txt"

    def trait2table(self, table_path):
        table_obj = Table()
        print "table_path", table_path
        table_obj.set_path(table_path)
        if table_obj.check():
            if table_obj.file_type == "xls":
                table_obj.get_xls_sheet1(table_path + '.xls')
                return table_path + '.xls'
            else:
                return table_path
        else:
            return False

    def check_trait(self, trait_path):
        with open(trait_path, 'rb') as f:
            def is_number(s):
                try:
                    float(s)
                    return True
                except ValueError:
                    pass

            '''
            header = f.readline()
            encoding_format = chardet.detect(header)
            if not (encoding_format['confidence'] == 1 and
                    encoding_format['encoding'].lower().startswith(('ascii', 'utf-8'))):
                return json.dumps({'success': False, 'info': "文件编码格式必须是ascii或utf-8", 'code': 'C2903402', 'variables': ''})
            '''
            phenoypes = dict()
            phenoype_list = list()
            # first_sample, test_one = f.readline().strip().split()[0:2]
            for phenoype in f.readline().strip().split()[1:]:
                phenoypes[phenoype] = list()
                phenoype_list.append(phenoype)
            if self.input_data.trait_type.lower() == "discrete":
                if len(phenoype_list) >= 2:
                    return json.dumps({'success': False, 'info': "选择非连续类型时，表型数据只可以为1列", "code": "C2903610"})
            trait_samples = list()
            for line in f:
                cols = line.strip().split()
                if len(cols) < len(phenoype_list) + 1:
                    return json.dumps(
                        {'success': False, 'info': "样本%s表型数据缺失", "variables": [cols[0]], "code": "C2903611"})
                for n, p in enumerate(phenoype_list):
                    phe = cols[n + 1]
                    if self.input_data.trait_type.lower() == "discrete":
                        if is_number(phe):
                            return json.dumps(
                                {'success': False, 'info': "选择非连续类型时，表型数据不应该是数字", 'code': 'C2903612', 'variables': ''})
                    else:
                        if not is_number(phe):
                            return json.dumps(
                                {'success': False, 'info': "选择连续类型时，表型数据必须是数字", 'code': 'C2903613', 'variables': ''})
                    phenoypes[p].append(phe)
                trait_samples.append(cols[0])
            select_samples = self.samples
            if len(trait_samples) != len(set(trait_samples)):
                return json.dumps({'success': False, 'info': "每个样本表型数据只可以有一行", "code": "C2903614"})
            if len(trait_samples) != len(select_samples):
                print trait_samples, select_samples
                return json.dumps({'success': False, 'info': "表型样本与预处理样本不一致！", 'code': 'C2903615', 'variables': ''})

            for p, phe in phenoypes.items():
                if len(set(phe)) == 1:
                    return json.dumps(
                        {'success': False, 'info': "表型%s只有单一的值请删除该列", "variables": [p], "code": "C2903616"})
        return True

    @check_sig
    def POST(self):
        data = web.input()
        self.input_data = data
        print data
        basic_args = ["task_id", "submit_location", 'task_type']
        basic_args += ['level', 'wgcna_module_id', 'trait_path', 'corr_method', 'trait_type', "file_id"]
        # check arg
        for arg in basic_args:
            if not hasattr(data, arg):
                variables = []
                variables.append(arg)
                info = {'success': False, 'info': "Lack argument: %s" % arg, 'code': 'C2903601', 'variables': variables}
                return json.dumps(info)
        main_info = self.medical_transcriptome.get_main_info(data.wgcna_module_id, 'sg_wgcna_module', data.task_id)
        project_sn = main_info["project_sn"]
        self.samples = main_info['samples']
        task_id = data.task_id
        '''
        # get abs path of input file and check
        client = data.client if hasattr(data,"client") else web.ctx.env.get('HTTP_CLIENT')
        if client == 'client01':
            target_dir = 'data'
        else:
            target_dir = 'tsanger-data'

        base_path = "/mnt/ilustre/{}/".format(target_dir)
        '''
        trait_path = data.trait_path
        if os.path.exists(trait_path):
            pass
        elif re.match(r'rerewrweset', data.trait_path):
            inter_dir = self.create_tmp_dir(data.task_id, "sg_wgcna_pipeline/")
            trait_path = self.download_from_s3(trait_path, inter_dir=inter_dir)
        elif re.match(r'^\w+://\S+/.+$', data.trait_path):
            inter_dir = self.create_tmp_dir(data.task_id, "sg_wgcna_pipeline/")
            trait_path = self.download_from_s3(data.trait_path, inter_dir=inter_dir)
        else:
            raise "文件传递格式错误 {}".format(data.trait_path)

        self.trait_path_clean = trait_path
        clean_path = self.trait2table(trait_path)
        print clean_path
        if clean_path:
            self.trait_path_clean = clean_path
        else:
            return json.dumps({'success': False, 'info': "不支持除了Unicode文本和excel之外的格式", "code": "C2903617"})
        check_trait = self.check_trait(self.trait_path_clean)
        if check_trait == True:
            pass
        else:
            return check_trait
        self.trait_path_clean = self.clean_trait(self.trait_path_clean)

        # create main table record
        params = dict(
            task_id=task_id,
            submit_location=data.submit_location,
            task_type=int(data.task_type),
            wgcna_module_id=data.wgcna_module_id,
            level=data.level,
            corr_method=data.corr_method,
            trait_type=data.trait_type,
            trait_path=data.trait_path,
            file_id=data.file_id,
        )
        params = json.dumps(params, sort_keys=True, separators=(',', ':'))

        if isinstance(data.wgcna_module_id, types.StringTypes):
            wgcna_module_id = ObjectId(data.wgcna_module_id)
        elif isinstance(data.wgcna_module_id, ObjectId):
            wgcna_module_id = data.wgcna_module_id
        else:
            raise Exception("wgcna_prepare_id参数必须为字符串或者ObjectId类型!")
        module_info = self.medical_transcriptome.get_main_info_by_record('sg_wgcna_module', main_id=wgcna_module_id,
                                                                       task_id=data.task_id)
        module_rdata = os.path.join(module_info["output_dir"], 'blockwiseModules_result.RData')

        name = "Relate" + '_'
        time_now = datetime.datetime.now()
        name += time_now.strftime("%Y%m%d_%H%M%S")

        main_info = dict(
            project_sn=project_sn,
            task_id=task_id,
            name=name,
            version="v1",
            created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
            wgcna_module_id=wgcna_module_id,
            desc='wgcna relation analysis main table',
            level=data.level,
            params=params,
            status="start"
        )
        main_id = self.medical_transcriptome.insert_main_table('sg_wgcna_relate', main_info)
        new_task_id = self.medical_transcriptome.get_new_id(data.task_id)
        main_table_data = {'run_id': new_task_id}

        # prepare option for workflow
        options = {
            "exp_eigengenes": data.wgcna_module_id,
            "main_id": str(main_id),
            "trait_path": self.trait_path_clean,
            "trait_type": data.trait_type,
            "corr_method": data.corr_method,
            "block_Rdata": module_rdata,
            "level": data.level,
            'main_table_data': main_table_data,
            "update_info": json.dumps({str(main_id): "sg_wgcna_relate"})  # to update status
        }
        # prepare to file
        to_files = ["medical_transcriptome.export_wgcna_relate_input(exp_eigengenes)"]

        # 把参数交给workflow运行相应的tool， 其中to_file用于准备tool的输入文件
        task_name = 'medical_transcriptome.report.wgcna_relate'
        self.set_sheet_data(name=task_name,
                            options=options,
                            main_table_name=name,  # 设置交互分析结果目录名
                            module_type="workflow",
                            to_file=to_files,
                            project_sn=project_sn,
                            new_task_id=new_task_id,
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
            _ = self.medical_transcriptome.update_group_is_use(data.task_id, data.group_id)
        if 'control_id' in data:
            _ = self.medical_transcriptome.update_group_compare_is_use(data.task_id, data.control_id)
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
        cmd += "s/medical_transcriptome/wgcna_relate "
        cmd += "-b http://bcl.tsg.com "
        args = dict(
            task_id="medical_transcriptome",
            task_type="2",
            submit_location="wgcna_relate",
            wgcna_module_id="5f583bf317b2bf35750c29d1",
            level="G",
            trait_path="/mnt/ilustre/users/sanger-dev/workspace/tmp/tmp_s3/tsg_38283/wgcna_pipeline/trait.txt.checked.txt",
            file_id="aaa",
            trait_type="discrete",
            corr_method="pearson",
        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)


if __name__ == '__main__':
    unittest.main()

# -*- coding: utf-8 -*-
import web
import json
import datetime
from collections import OrderedDict
from mainapp.controllers.project.ref_rna_controller import RefRnaController
from mbio.api.to_file.ref_rna import *
from mainapp.libs.signature import check_sig
import unittest
import os,re
import chardet


class WgcnaPipelineAction(RefRnaController):
    def __init__(self):
        super(WgcnaPipelineAction, self).__init__(instant=False)
        self.task_name = 'ref_rna.report.wgcna_pipeline'
        self.expected_args = ["task_id", "submit_location", 'task_type']
        self.expected_args += ['group_id', 'group_dict', 'exp_id', 'me', 'geneset_id', 'exp_level', 'cv']
        self.expected_args += ['mergeCutHeight', 'power', 'minModuleSize', 'networkType', 'minKMEtoStay']
        self.expected_args += ['trait_path', 'corr_method', 'trait_type', "file_id"]
        self.expected_args += ['threshold', 'top']
        self.input_data = web.input()
        self.collection = "sg_wgcna_pipeline"

    @check_sig
    def POST(self):
        check = self.check_params()
        if check is not True:
            return check
        packed_params = self.pack_params()
        main_id, main_table_name = self.create_main_table(packed_params, self.collection)
        self.prepare_workflow_options(main_id, self.collection, main_table_name)
        task_info = super(WgcnaPipelineAction,self).POST()
        task_info['content'] = {'ids': {'id': str(main_id), 'name': main_table_name}}
        return json.dumps(task_info)

    def check_params(self):
        for arg in self.expected_args:
            if not hasattr(self.input_data, arg):
                info = {'success': False, 'info': "Lack argument: {}".format(arg)}
                return json.dumps(info)
        # do some other check
        target_dir = 'data' if self.input_data.client == "client01" else 'tsanger-data'
        base_path = "/mnt/ilustre/{}/".format(target_dir)
        trait_path = base_path + self.input_data.trait_path
        if os.path.exists(trait_path):
            pass
        elif re.match(r'^\w+://\S+/.+$', self.input_data.trait_path) or re.match(r'/mnt/ilustre', self.input_data.trait_path):
            inter_dir = self.create_tmp_dir(self.input_data.task_id, "wgcna_pipeline/")
            trait_path = self.download_from_s3(self.input_data.trait_path, inter_dir=inter_dir)
        elif self.input_data.trait_path.startswith("rerewrweset"):
            inter_dir = self.create_tmp_dir(self.input_data.task_id, "wgcna_pipeline/")
            trait_path = self.download_from_s3("/mnt/ilustre/data/" + self.input_data.trait_path, inter_dir=inter_dir)
        else:
            raise Exception("文件传递格式错误 {}".format(self.input_data.trait_path))

        self.trait_path = trait_path
        with open(trait_path, 'rb') as f:
            def is_number(s):
                try:
                    float(s)
                    return True
                except ValueError:
                    pass
            header = f.readline()
            encoding_format = chardet.detect(header)
            if not (encoding_format['confidence'] == 1 and
                    encoding_format['encoding'].lower().startswith(('ascii', 'utf-8'))):
                return json.dumps({'success': False, 'info': "文件编码格式必须是ascii或utf-8"})
            first_sample, test_one = f.readline().strip().split()[0:2]
            if self.input_data.trait_type.lower() == "discrete":
                if is_number(test_one):
                    return json.dumps({'success': False, 'info': "选择非连续类型时，表型数据不应该是数字"})
            else:
                if not is_number(test_one):
                    return json.dumps({'success': False, 'info': "选择连续类型时，表型数据必须是数字"})
            trait_samples = [first_sample] + [x.strip().split()[0] for x in f if x.strip()]
            select_samples = [y for x in json.loads(self.input_data.group_dict).values() for y in x]
            if len(trait_samples) != len(select_samples):
                return json.dumps({'success': False, 'info': "表型样本与预处理样本不一致！"})
        return True

    def pack_params(self):
        input_data_dict = dict(self.input_data)
        params_dict = dict()
        for each in self.expected_args:
            if each == "task_type":
                params_dict[each] = int(input_data_dict[each])
            elif each == "group_dict":
                params_dict[each] = json.loads(input_data_dict[each], object_pairs_hook=OrderedDict)
            else:
                params_dict[each] = input_data_dict[each]
        return json.dumps(params_dict, sort_keys=True, separators=(',',':'))

    def create_main_table(self, packed_params, collection_name):
        """
        对应的workflow会调用4个tool, 4个tool对应四次交互分析结果，即对应四个主表，他们用pipeline_id追踪来自哪一次的一键化分析。
        """
        result_info = self.ref_rna.get_main_info(self.input_data.exp_id, 'sg_express')
        self.project_sn = result_info["project_sn"]
        exp_info = json.loads(result_info['params'])
        exp_level = self.input_data.exp_level[0].upper()
        exp_type = exp_info["type"]
        quant_method = exp_info["express_method"]
        name = "WgcnaPipeline" + '_' + exp_level + '_' + quant_method + '_' + exp_type.upper() + '_'
        time_now = datetime.datetime.now()
        name += time_now.strftime("%Y%m%d_%H%M%S")
        main_info = dict(
            project_sn=self.project_sn,
            task_id=self.input_data.task_id,
            name=name,
            created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
            exp_id=self.input_data.exp_id,
            desc='wgcna pipeline analysis main table',
            type=self.input_data.exp_level,
            params=packed_params,
            status="start")
        main_id = self.ref_rna.insert_main_table(collection_name, main_info)
        if self.input_data.geneset_id.lower() not in ["all", "refall", "none"]:
            self.ref_rna.insert_geneset_info(self.input_data.geneset_id, collection_name, str(main_id))
        return main_id, name

    def prepare_workflow_options(self, main_id, collection_name, main_table_name):
        input_data_dict = dict(self.input_data)
        options = dict()
        for each in self.expected_args:
            options[each] = input_data_dict[each]
        group_dict = json.loads(self.input_data.group_dict, object_pairs_hook=OrderedDict)
        if str(self.input_data.group_id).lower() == 'all':
            samples = group_dict['all']
            group_dict = OrderedDict([(x, [x]) for x in samples])
        options.update({
            "trait_path": self.trait_path,
            "main_id": str(main_id),
            "exp_matrix": self.input_data.exp_id+";"+self.input_data.geneset_id,
            "group_dict": json.dumps(group_dict),
            "raw_group_dict": self.input_data.group_dict,
            "update_info": json.dumps({str(main_id): collection_name})  # to update sg_status
        })
        # prepare to file
        to_files = ["ref_rna.export_geneset_exp_matrix(exp_matrix)"]
        # 把参数交给workflow运行相应的tools， 其中to_file用于准备tool的输入文件
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
        cmd += "s/ref_rna/wgcna_pipeline "
        cmd += "-b http://192.168.12.102:9090 "
        args = dict(
            task_id="tsg_28226",
            task_type="2",
            submit_location="WgcnaPipeline",
            exp_id="5a9b68d0a4e1af16f6d8e47b",
            geneset_id="all",
            group_id="5a9b685aa4e1af16f6d60f7c",
            group_dict=r'{"A": ["5a9b6857a4e1af16f6d5fa56", "5a9b6857a4e1af16f6d5fa57", "5a9b6857a4e1af16f6d5fa58"], "B": ["5a9b6857a4e1af16f6d5fa59", "5a9b6857a4e1af16f6d5fa5a", "5a9b6857a4e1af16f6d5fa5b"], "C": ["5a9b6857a4e1af16f6d5fa5c", "5a9b6857a4e1af16f6d5fa5d", "5a9b6857a4e1af16f6d5fa5e"], "D": ["5a9b6857a4e1af16f6d5fa60", "5a9b6857a4e1af16f6d5fa61", "5a9b6857a4e1af16f6d5fa63"], "E": ["5a9b6857a4e1af16f6d5fa55", "5a9b6857a4e1af16f6d5fa5f", "5a9b6857a4e1af16f6d5fa62"], "F": ["5a9b6857a4e1af16f6d5fa52", "5a9b6857a4e1af16f6d5fa53", "5a9b6857a4e1af16f6d5fa54"]}'.replace('"', '\\"'),
            me='2.10',
            cv='0.6',
            exp_level="gene",
            mergeCutHeight="0.45",
            power="7",
            minModuleSize="30",
            networkType="signed",
            minKMEtoStay="0.5",
            trait_path="/rerewrweset/files/m_188/188_59f029eed6609/tsg_28226/trait_data/traits_1522310047.txt",
            trait_type="discrete",
            corr_method="pearson",
            file_id="6919032",
            threshold="0.03",
            top="30",
        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)


if __name__ == '__main__':
    unittest.main()

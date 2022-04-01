# -*- coding: utf-8 -*-

import datetime
import json
import os
import re
import unittest
from collections import OrderedDict
import web
from mainapp.controllers.project.medical_transcriptome_controller import MedicalTranscriptomeController
from mainapp.libs.signature import check_sig
from mbio.packages.rna.table_check import Table
from mbio.api.to_file.medical_transcriptome import *
from bson.objectid import ObjectId


class WgcnaPipelineAction(MedicalTranscriptomeController):
    def __init__(self):
        super(WgcnaPipelineAction, self).__init__(instant=False)
        self.task_name = 'medical_transcriptome.report.wgcna_pipeline'
        self.expected_args = ["task_id", "submit_location", 'task_type']
        self.expected_args += ['exp_id', 'group_id', 'group_dict', 'me', 'geneset_id', 'level', 'cv']
        self.expected_args += ['mergeCutHeight', 'power', 'minModuleSize', 'networkType', 'minKMEtoStay']
        self.expected_args += ['trait_path', 'corr_method', 'trait_type', "file_id"]
        self.expected_args += ['threshold', 'top']
        self.input_data = web.input()
        self.collection = "sg_wgcna_pipeline"

    @check_sig
    def POST(self):
        print self.input_data
        check = self.check_params()
        if check is not True:
            return check
        packed_params = self.pack_params()
        main_id, main_table_name = self.create_main_table(packed_params, self.collection)
        self.prepare_workflow_options(main_id, self.collection, main_table_name)
        task_info = super(WgcnaPipelineAction, self).POST()
        task_info['content'] = {'ids': {'id': str(main_id), 'name': main_table_name}}
        if 'group_id' in self.input_data and str(self.input_data.group_id).lower() != 'all':
            _ = self.medical_transcriptome.update_group_is_use(self.input_data.task_id, self.input_data.group_id)
            '''
            group_id 根据类型联动
            if self.input_data.category in ['mRNA', 'lncRNA']:
                _ = self.whole_transcriptome.update_group_is_use(self.input_data.task_id, self.input_data.group_id)
            elif self.input_data.category  == "miRNA":
                _ = self.whole_transcriptome.update_group_is_use(self.input_data.task_id, self.input_data.group_id)
            elif self.input_data.category  == "circRNA":
                _ = self.whole_transcriptome.update_group_is_use(self.input_data.task_id, self.input_data.group_id)
            '''

        if 'control_id' in self.input_data:
            _ = self.medical_transcriptome.update_group_compare_is_use(self.input_data.task_id,
                                                                     self.input_data.control_id)
        return json.dumps(task_info)

    '''
    NUM_LETTER = re.compile("^(?!\d+$)[\da-zA-Z_]+$")     #数字和字母组合，不允许纯数字
    FIRST_LETTER = re.compile("^[a-zA-Z]")           #只能以字母开头

    def account_name_fomart(Name):
        if NUM_LETTER.search(Name):
           if FIRST_LETTER.search(Name):
               return True
        return False

    def is_account_name_fomart(Name):
            if not account_name_fomart(Name):
            msg = "用户名不合法"
            return msg
    '''

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
                    return json.dumps({'success': False, 'info': "选择非连续类型时，表型数据只可以为1列", "code": "C2903411"})
            trait_samples = list()
            for line in f:
                cols = line.strip().split()
                if len(cols) < len(phenoype_list) + 1:
                    return json.dumps(
                        {'success': False, 'info': "样本%s表型数据缺失", "variables": [cols[0]], "code": "C2903412"})
                for n, p in enumerate(phenoype_list):
                    phe = cols[n + 1]
                    if self.input_data.trait_type.lower() == "discrete":
                        if is_number(phe):
                            return json.dumps(
                                {'success': False, 'info': "选择非连续类型时，表型数据不应该是数字", 'code': 'C2903403', 'variables': ''})
                    else:
                        if not is_number(phe):
                            return json.dumps(
                                {'success': False, 'info': "选择连续类型时，表型数据必须是数字", 'code': 'C2903404', 'variables': ''})
                    phenoypes[p].append(phe)
                trait_samples.append(cols[0])
            select_samples = [y for x in json.loads(self.input_data.group_dict).values() for y in x]
            if len(trait_samples) != len(set(trait_samples)):
                return json.dumps({'success': False, 'info': "每个样本表型数据只可以有一行", "code": "C2903413"})
            if len(trait_samples) != len(select_samples):
                print trait_samples, select_samples
                return json.dumps({'success': False, 'info': "表型样本与预处理样本不一致！", 'code': 'C2903405', 'variables': ''})

            for p, phe in phenoypes.items():
                if len(set(phe)) == 1:
                    return json.dumps(
                        {'success': False, 'info': "表型%s只有单一的值请删除该列", "variables": [p], "code": "C2903414"})
        return True

    def clean_trait(self, trait_path):
        with open(trait_path, 'rb') as f, open(trait_path + ".checked.txt", 'w') as fo:
            header = f.readline()
            cols = header.strip().lstrip("#").split()
            fo.write("\t".join(cols) + "\n")
            for line in f:
                cols = line.strip().split()
                fo.write("\t".join(cols) + "\n")
        return trait_path + ".checked.txt"

    def check_params(self):
        for arg in self.expected_args:
            if not hasattr(self.input_data, arg):
                variables = []
                variables.append(arg)
                info = {'success': False, 'info': "Lack argument: %s" % arg, 'code': 'C2903401', 'variables': variables}
                return json.dumps(info)
        # do some other check
        target_dir = 'data' if self.input_data.client == "client01" else 'tsanger-data'
        # base_path = "/mnt/ilustre/{}/".format(target_dir)

        trait_path = self.input_data.trait_path
        if os.path.exists(trait_path):
            pass
        elif re.match(r'^\w+://\S+/.+$', self.input_data.trait_path) or re.match(r'/mnt/ilustre',
                                                                                 self.input_data.trait_path):
            inter_dir = self.create_tmp_dir(self.input_data.task_id, "wgcna_pipeline/")
            trait_path = self.download_from_s3(self.input_data.trait_path, inter_dir=inter_dir)
        elif self.input_data.trait_path.startswith("rerewrweset"):
            inter_dir = self.create_tmp_dir(self.input_data.task_id, "wgcna_pipeline/")
            trait_path = self.download_from_s3("/mnt/ilustre/data/" + self.input_data.trait_path, inter_dir=inter_dir)
        elif re.match(r'rerewrweset', self.input_data.trait_path):
            inter_dir = self.create_tmp_dir(self.input_data.task_id, "wgcna_pipeline/")
            trait_path = self.download_from_s3(self.input_data.trait_path, inter_dir=inter_dir)
        else:
            raise Exception("文件传递格式错误 {}".format(self.input_data.trait_path))
        self.trait_path_clean = trait_path
        clean_path = self.trait2table(trait_path)
        print clean_path
        if clean_path:
            self.trait_path_clean = clean_path
        else:
            return json.dumps({'success': False, 'info': "不支持除了Unicode文本和excel之外的格式", "code": "C2903415"})

        check_trait = self.check_trait(self.trait_path_clean)
        if check_trait == True:
            pass
        else:
            return check_trait

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
        return json.dumps(params_dict, sort_keys=True, separators=(',', ':'))

    def create_main_table(self, packed_params, collection_name):
        """
        对应的workflow会调用4个tool, 4个tool对应四次交互分析结果，即对应四个主表，他们用pipeline_id追踪来自哪一次的一键化分析。
        """

        self.exp_info = self.medical_transcriptome.get_main_info(self.input_data.exp_id, 'sg_exp', self.input_data.task_id)
        self.project_sn = self.exp_info["project_sn"]
        '''
        exp_info = json.loads(self.exp_info['params'])
        exp_level = self.input_data.exp_level[0].upper()
        exp_type = exp_info["level"]
        '''
        exp_level = self.exp_info["level"]
        exp_type = self.exp_info["exp_type"]
        try:
            quant_method = self.exp_info["method"]
        except:
            quant_method = "tpm"
        # exp_info["method"]
        # name = "WgcnaPipeline" + '_' + self.input_data.category + '_' + quant_method + '_' + exp_type.upper() + '_'
        name = "Pipeline" + '_' + exp_level + '_' + quant_method + '_' + exp_type + '_'
        time_now = datetime.datetime.now()
        name += time_now.strftime("%Y%m%d_%H%M%S")
        main_info = dict(
            project_sn=self.project_sn,
            task_id=self.input_data.task_id,
            name=name,
            version="v1",
            created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
            exp_id=str(self.exp_info['main_id']),
            desc='wgcna pipeline analysis main table',
            level=self.input_data.level,
            # category=self.input_data.category,
            params=packed_params,
            status="start")
        main_id = self.medical_transcriptome.insert_main_table(collection_name, main_info)
        if self.input_data.geneset_id.lower() not in ["all", "refall", "none"]:
            self.medical_transcriptome.insert_geneset_info(self.input_data.geneset_id, collection_name, str(main_id))
        return main_id, name

    def prepare_workflow_options(self, main_id, collection_name, main_table_name):
        self.trait_path_clean = self.clean_trait(self.trait_path_clean)
        input_data_dict = dict(self.input_data)
        options = dict()
        for each in self.expected_args:
            options[each] = input_data_dict[each]
        group_dict = json.loads(self.input_data.group_dict, object_pairs_hook=OrderedDict)
        if str(self.input_data.group_id).lower() == 'all':
            samples = group_dict['all']
            group_dict = OrderedDict([(x, [x]) for x in samples])
        connect_exp = self.db['sg_exp']
        record_exp = connect_exp.find_one({'task_id': self.input_data.task_id, 'main_id': ObjectId(self.input_data.exp_id)})
        is_rmbe = str(record_exp['is_rmbe']).lower()
        if is_rmbe == 'false':
            exp_id = self.input_data.exp_id
        if is_rmbe == 'true':
            exp_id = str(record_exp['batch_main_id'])
        new_task_id = self.medical_transcriptome.get_new_id(self.input_data.task_id)
        main_table_data = {'run_id': new_task_id}
        options.update({
            "trait_path": self.trait_path_clean,
            "main_id": str(main_id),
            # "exp_id": str(self.exp_info['main_id']),
            "exp_matrix": str(exp_id) + ';' + self.input_data.geneset_id + ';' +
                          self.input_data.level + ';' + is_rmbe,
            "gene_name": str(exp_id) + ";" + self.input_data.geneset_id,
            "group_dict": json.dumps(group_dict),
            "raw_group_dict": self.input_data.group_dict,
            "main_table_data": main_table_data,
            "update_info": json.dumps({str(main_id): collection_name})  # to update sg_status
        })
        # prepare to file
        to_files = ["medical_transcriptome.export_geneset_exp_matrix_new(exp_matrix)",
                    "medical_transcriptome.export_geneset_genes_name(gene_name)"
                    ]

        # 把参数交给workflow运行相应的tools， 其中to_file用于准备tool的输入文件
        self.set_sheet_data(
            name=self.task_name,
            options=options,
            main_table_name=main_table_name,  # 设置交互分析结果目录名
            module_type="workflow",
            to_file=to_files,
            project_sn=self.project_sn,
            new_task_id=new_task_id,
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
        cmd += "s/medical_transcriptome/wgcna_pipeline "
        cmd += "-b http://bcl.tsg.com "
        args = dict(
            task_id="medical_transcriptome",
            task_type="2",
            exp_id="5f50cacf17b2bf5a6c8bfd88",
            group_id='5f46228c17b2bf20e4e269e1',
            level='G',
            group_dict=json.dumps({"H1": ["H1581_1", "H1581_2", 'H1581_3'], "H2": ["H1581_4", "H1581_5", 'H1581_6'], "H3": ["H1581_7", "H1581_8", 'H1581_9'], 'S1':["SNU16_1", "SNU16_2", "SNU16_3"], 'S2':["SNU16_4", "SNU16_5", "SNU16_6"], 'S3':["SNU16_7", "SNU16_8", "SNU16_9"]}).replace('"', '\\"'),
            submit_location="wgcnapipeline",
            geneset_id="All",
            me='1',
            cv='0.1',
            mergeCutHeight="0.45",
            power="0",
            minModuleSize="30",
            networkType="signed",
            minKMEtoStay="0.3",
            trait_path="/mnt/ilustre/users/sanger-dev/workspace/tmp/tmp_s3/tsg_38283/wgcna_pipeline/trait.txt.checked.txt",
            trait_type="discrete",
            corr_method="pearson",
            file_id="6919032",
            threshold="0.2",
            top="30",
        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)


if __name__ == '__main__':
    unittest.main()

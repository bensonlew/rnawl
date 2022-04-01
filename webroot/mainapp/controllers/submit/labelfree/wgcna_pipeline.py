# -*- coding: utf-8 -*-
import web
import json
import datetime
from collections import OrderedDict
from mainapp.controllers.project.labelfree_controller import LabelfreeController
from mbio.api.to_file.labelfree import *
from mbio.packages.rna.table_check import Table
from mainapp.libs.signature import check_sig
import unittest
import os,re
import chardet


class WgcnaPipelineAction(LabelfreeController):
    def __init__(self):
        super(WgcnaPipelineAction, self).__init__(instant=False)
        self.task_name = 'labelfree.report.wgcna_pipeline'
        self.expected_args = ["task_id", "submit_location", 'task_type']
        self.expected_args += ['group_id', 'group_dict', 'exp_id', 'me', 'proteinset_id', 'cv']
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
        if 'group_id' in self.input_data and str(self.input_data.group_id).lower() != 'all':
            _ = self.labelfree.update_group_is_use(self.input_data.task_id, self.input_data.group_id)
        if 'control_id' in self.input_data:
            _ = self.labelfree.update_group_compare_is_use(self.input_data.task_id, self.input_data.control_id)
        return json.dumps(task_info)

    def trait2table(self, table_path):
        table_obj = Table()
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
                    return json.dumps({'success': False, 'info': "选择非连续类型时，表型数据只可以为1列"})
            trait_samples = list()
            for line in f:
                cols = line.strip().split()
                if len(cols) < len(phenoype_list) + 1:
                    return json.dumps({'success': False, 'info': "样本{}表型数据缺失".format(';'.join(cols))})
                for n, p in enumerate(phenoype_list):
                    phe = cols[n+1]
                    if self.input_data.trait_type.lower() == "discrete":
                        if is_number(phe):
                            return json.dumps({'success': False, 'info': "选择非连续类型时，表型数据不应该是数字", 'code': 'C2903403', 'variables': ''})
                    else:
                        if not is_number(phe):
                            return json.dumps({'success': False, 'info': "选择连续类型时，表型数据必须是数字", 'code': 'C2903404', 'variables': ''})
                    phenoypes[p].append(phe)
                trait_samples.append(cols[0])
            select_samples = [y for x in json.loads(self.input_data.group_dict).values() for y in x]
            if len(trait_samples) != len(set(trait_samples)):
                return json.dumps({'success': False, 'info': "每个样本表型数据只可以有一行"})
            if len(trait_samples) != len(select_samples):
                return json.dumps({'success': False, 'info': "表型样本与预处理样本不一致！", 'code': 'C2903405', 'variables': ''})

            for p, phe in phenoypes.items():
                if len(set(phe)) == 1:
                    return json.dumps({'success': False, 'info': "表型{}只有单一的值请删除该列".format(p)})
        return True
    def clean_trait(self, trait_path):
        with open(trait_path, 'rb') as f, open(trait_path+ ".checked.txt", 'w') as fo:
            header = f.readline()
            cols = header.strip().lstrip("#").split()
            fo.write("\t".join(cols) + "\n")
            for line in f:
                cols = line.strip().split()
                fo.write("\t".join(cols) + "\n")
        return trait_path+ ".checked.txt"

    def check_params(self):
        for arg in self.expected_args:
            if not hasattr(self.input_data, arg):
                variables = []
                variables.append(arg)
                info = {'success': False, 'info': "Lack argument: %s" % arg, 'code': 'C2903401', 'variables': variables}
                return json.dumps(info)
        # do some other check
        target_dir = 'data' if self.input_data.client == "client01" else 'tsanger-data'
        base_path = "/mnt/ilustre/{}/".format(target_dir)

        trait_path = base_path + self.input_data.trait_path
        # trait_path = self.input_data.trait_path
        if os.path.exists(trait_path):
            pass
        elif re.match(r'^\w+://\S+/.+$', self.input_data.trait_path) or re.match(r'/mnt/ilustre', self.input_data.trait_path):
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
        if clean_path:
            self.trait_path_clean = clean_path
        else:
            return json.dumps({'success': False, 'info': "不支持除了Unicode文本和excel之外的格式"})

        os.system("sed -i -e '/^\s*$/d' -e 's/\xEF\xBB\xBF//' {}".format(self.trait_path_clean))#删除空行和bom信息
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
        return json.dumps(params_dict, sort_keys=True, separators=(',',':'))

    def create_main_table(self, packed_params, collection_name):
        """
        对应的workflow会调用4个tool, 4个tool对应四次交互分析结果，即对应四个主表，他们用pipeline_id追踪来自哪一次的一键化分析。
        """

        result_info = self.labelfree.get_main_info(self.input_data.exp_id, 'sg_express', self.input_data.task_id)
        self.project_sn = result_info["project_sn"]
        name = "WgcnaPipeline" + '_'
        time_now = datetime.datetime.now()
        name += time_now.strftime("%Y%m%d_%H%M%S")
        main_info = dict(
            project_sn=self.project_sn,
            task_id=self.input_data.task_id,
            name=name,
            version="v2",
            created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
            exp_id=self.input_data.exp_id,
            desc='wgcna pipeline analysis main table',
            params=packed_params,
            status="start")
        main_id = self.labelfree.insert_main_table(collection_name, main_info)
        if self.input_data.proteinset_id.lower() not in ["all", "none"]:
            self.labelfree.insert_proteinset_info(self.input_data.proteinset_id, collection_name, str(main_id))
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
        options.update({
            "trait_path": self.trait_path_clean,
            "main_id": str(main_id),
            "exp_matrix": self.input_data.exp_id+";"+self.input_data.proteinset_id,
            "group_dict": json.dumps(group_dict),
            "raw_group_dict": self.input_data.group_dict,
            "update_info": json.dumps({str(main_id): collection_name})  # to update sg_status
        })
        # prepare to file
        to_files = ["labelfree.export_proteinset_exp_matrix(exp_matrix)"]
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
        cmd += "s/labelfree/wgcna_pipeline "
        cmd += "-b http://bcl.tsg.com "
        args = dict(
            task_id="tsg_34739",
            task_type="2",
            exp_id='5d22a82817b2bf1c689dbed6',
            group_id='5d22a82717b2bf1c689db738',
            group_dict=r'{"F":["F_1","F_2","F_3"],"S":["S_1","S_2","S_3"]}'.replace(
                '"', '\\"'),
            submit_location="wgcnapipeline",
            proteinset_id="all",
            me='1.50',
            cv='0.1',
            mergeCutHeight="0.45",
            power="7",
            minModuleSize="30",
            networkType="signed",
            minKMEtoStay="0.5",
            trait_path="/mnt/ilustre/users/sanger-dev/sg-users/deqing/projects_test_input/traits_fyt.xls",
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

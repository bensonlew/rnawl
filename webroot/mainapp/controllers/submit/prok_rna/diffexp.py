# -*- coding: utf-8 -*-
import web
import json
import datetime
from collections import OrderedDict
from mainapp.controllers.project.prok_rna_controller import ProkRNAController
from mbio.api.to_file.prok_rna import *
from mainapp.libs.signature import check_sig
import os
import unittest


class DiffexpAction(ProkRNAController):
    def __init__(self):
        super(DiffexpAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        basic_args = ["task_id", "submit_location", 'task_type']
        basic_args += ['group_id', 'group_dict', 'control_id',
                       'exp_id', 'exp_level', 'diff_method',
                       'stat_type', 'stat_cutoff', 'fc']
        # check arg
        for arg in basic_args:
            if not hasattr(data, arg):
                info = {'success': False, 'info': "Lack argument: {}".format(arg)}
                return json.dumps(info)
        exp_info = self.prok_rna.get_exp_params_info(data.exp_id, data.task_id)
        project_sn = exp_info["project_sn"]
        task_id = data.task_id
        group_dict = json.loads(data.group_dict, object_pairs_hook=OrderedDict)
        if str(data.group_id).lower() == 'all':
            samples = group_dict['all']
            group_dict = OrderedDict([(x, [x]) for x in samples])
        group_size = list()
        success = list()
        for key in group_dict:
            group_size.append(len(group_dict[key]))
        group_size.sort()
        if group_size[0] == group_size[-1] == 1:
            if data.diff_method == "DESeq2":
                success.append('DESeq2不适合处理单样本和单样本比较， 请选择DEGseq或edgeR')
        # elif group_size[0] == 1 and group_size[1] >= 2:
        #     if data.diff_method == "DESeq2" or data.diff_method == "DEGseq":
        #         success.append('你的分组方案表明可能要进行单样本和多样本的比较，'
        #                        '此时请选择edgeR做差异分析或者重新设计分组方案')
        elif group_size[0] >= 2:
            if data.diff_method == "DEGseq":
                success.append('只涉及包含多样本的组与组间比较时，我们只推荐DESeq2或者edgeR')
        if success:
            info = {'success': False, 'info': success[0]}
            return json.dumps(info)

        # create main table record
        exp_level = data.exp_level
        exp_type, quant_method = exp_info['exp_type'], exp_info['method']
        params = dict(
            task_id=task_id,
            submit_location=data.submit_location,
            task_type=int(data.task_type),
            exp_id=data.exp_id,
            group_id=data.group_id,
            control_id=data.control_id,
            exp_level=exp_level,
            group_dict=json.loads(data.group_dict, object_pairs_hook=OrderedDict),
            fc=data.fc,
            stat_type=data.stat_type,
            stat_cutoff=data.stat_cutoff,
            # quant_method=quant_method,
            diff_method=data.diff_method,
            # type=data.type,
        )
        # for special args
        if hasattr(data, 'correct_method'):
            correct_method = data.correct_method
            params.update(dict(correct_method=correct_method))
        else:
            correct_method = "BH"

        params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        name = "Diff" + '_' + exp_level + '_' + quant_method + '_' + data.diff_method + '_'
        time_now = datetime.datetime.now()
        name += time_now.strftime("%Y%m%d_%H%M%S")
        if type(params) == dict:
            params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        main_info = dict(
            project_sn=project_sn,
            task_id=task_id,
            name=name,
            created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
            desc='differential analysis main table',
            exp_level=exp_level,
            group_dict=group_dict,
            group_id=data.group_id,
            exp_id=data.exp_id,
            params=params,
            status="start",
        )
        main_id = self.prok_rna.insert_main_table('sg_diff', main_info)

        # prepare option for workflow
        options = {
            "exp_matrix": data.exp_id,
            "group_dict": json.dumps(group_dict),
            "diff_main_id": str(main_id),
            "method": data.diff_method,
            "count": self.prok_rna.get_count_path(data.exp_id),
            "fc": float(data.fc),
            "group": json.dumps(group_dict),
            "cmp": data.control_id,
            "pvalue": float(data.stat_cutoff),
            "pvalue_padjust": data.stat_type,
            "padjust_way": correct_method,
            "update_info": json.dumps({str(main_id): "sg_diff"}),
            # "type": data.type,
            "task_id": data.task_id,
            "exp_level": exp_level,
        }

        # prepare to file
        to_files = ["prok_rna.export_exp_matrix_prok(exp_matrix)",
                    "prok_rna.export_group(group)",
                    "prok_rna.export_compare(cmp)"]

        # 把参数交给workflow运行相应的tool， 其中to_file用于准备tool的输入文件
        task_name = 'prok_rna.report.diffexp'
        self.set_sheet_data(name=task_name,
                            options=options,
                            main_table_name=name,  # 设置交互分析结果目录名
                            module_type="workflow",
                            to_file=to_files,
                            project_sn=project_sn,
                            task_id=data.task_id)

        # 运行workflow 并传回参数
        task_info = super(DiffexpAction, self).POST()
        task_info['content'] = {
            'ids': {
                'id': str(main_id),
                'name': name
                }
        }
        if 'group_id' in data and str(data.group_id).lower() != 'all':
            _ = self.prok_rna.update_group_is_use(data.task_id, data.group_id)
        if 'control_id' in data:
            _ = self.prok_rna.update_group_compare_is_use(data.task_id, data.control_id)
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
        cmd += "s/prok_rna/diffexp "
        cmd += "-b http://192.168.12.102:9090 "
        args = dict(
            task_id="prok_rna_srna",
            task_type="2",
            submit_location="diffexp",
            exp_id="5b877be2a4e1af3eb5578ab9",
            group_id="5b7269fe77b3f3b1133b4a71",
            exp_level='mRNA+sRNA',
            group_dict=r'{"WT": ["WT_1", "WT_2", "WT_3"], "rcsBKO": ["rcsBKO_1", "rcsBKO_2", "rcsBKO_3"]}'.replace('"', '\\"'),
            # type='ref',
            control_id="5b726abe77b3f3b1133b6dab",
            diff_method='DEGseq',
            stat_type='padjust',
            stat_cutoff='0.05',
            fc='1.5',
        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(str(x) for x in arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTest(TestFunction('test_this'))
    unittest.TextTestRunner(verbosity=2).run(suite)

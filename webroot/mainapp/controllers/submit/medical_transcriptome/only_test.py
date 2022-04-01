# -*- coding: utf-8 -*-
import web
import json
import datetime
from collections import OrderedDict
from mainapp.controllers.project.medical_transcriptome_controller import MedicalTranscriptomeController
from mbio.api.to_file.medical_transcriptome import *
from mainapp.libs.signature import check_sig
import os
import unittest


class OnlyTestAction(MedicalTranscriptomeController):
    def __init__(self):
        super(OnlyTestAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        print data
        basic_args = ["task_id", "submit_location", 'task_type']
        basic_args += ['group_id', 'group_dict', 'control_id',
                       'exp_id', 'exp_level', 'diff_method',
                       'stat_type', 'stat_cutoff', 'fc', 'type']
        # check arg
        for arg in basic_args:
            if not hasattr(data, arg):
                variables = []
                variables.append(arg)
                info = {'success': False, 'info': "Lack argument:   %s" % arg, 'code': 'C2900301', 'variables': variables}
                return json.dumps(info)
        if not hasattr(data, 'tpm_filter_threshold'):
            data.tpm_filter_threshold = "0"
        elif not data.tpm_filter_threshold:
            data.tpm_filter_threshold = "0"
        exp_info = self.medical_transcriptome.get_exp_params_info(data.exp_id, data.task_id)
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
                success.append('If your experiment has no biological replicates, please select DEGseq or edgeR')
        elif group_size[0] == 1 and group_size[1] >= 2:
            if data.diff_method == "DESeq2":
                success.append('If some groups have biological replicates and some groups no, please select edgeR')
        elif group_size[0] >= 2:
            if data.diff_method == "DEGseq":
                success.append('If your experiment has biological replicates, please select  DESeq2 or edgeR.')
        if success:
            variables = []
            variables.append(success[0])
            info = {'success': False, 'info': '%s' % success[0], 'code': 'C2900302', 'variables': variables}
            return json.dumps(info)

        # create main table record
        exp_level = exp_info['exp_level']
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
            fc=str(float(data.fc)),
            stat_type=data.stat_type,
            stat_cutoff=data.stat_cutoff,
            # quant_method=quant_method,
            diff_method=data.diff_method,
            type=data.type,
            tpm_filter_threshold=str(float(data.tpm_filter_threshold))
        )
        # for special args
        if hasattr(data, 'correct_method'):
            correct_method = data.correct_method
            params.update(dict(correct_method=correct_method))
        else:
            correct_method = "BH"

        params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        name = "DiffExpress" + '_' + exp_level + '_' + quant_method + '_' + data.diff_method + '_'
        time_now = datetime.datetime.now()
        name += time_now.strftime("%Y%m%d_%H%M%S")
        if type(params) == dict:
            params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        main_info = dict(
            project_sn=project_sn,
            task_id=task_id,
            version="v3",
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
        main_id = self.medical_transcriptome.insert_main_table('sg_diff', main_info)

        # prepare option for workflow
        options = {
            "exp_matrix": data.exp_id,
            "group_dict": json.dumps(group_dict),
            "diff_main_id": str(main_id),
            "method": data.diff_method,
            "count": self.medical_transcriptome.get_count_path(data.exp_id),
            "fc": float(data.fc),
            "group": json.dumps(group_dict),
            "cmp": data.control_id,
            "pvalue": float(data.stat_cutoff),
            "pvalue_padjust": data.stat_type,
            "padjust_way": correct_method,
            "update_info": json.dumps({str(main_id): "sg_diff"}),
            "type": data.type,
            "task_id": data.task_id,
            "exp_level": exp_level,
            "exp_type": exp_type,
            "tpm_filter_threshold":data.tpm_filter_threshold
        }

        # prepare to file
        to_files = ["medical_transcriptome.export_exp_matrix(exp_matrix)",
                    "medical_transcriptome.export_group(group)",
                    "medical_transcriptome.export_compare(cmp)"]

        # 把参数交给workflow运行相应的tool， 其中to_file用于准备tool的输入文件
        task_name = 'medical_transcriptome.report.diffexp'
        self.set_sheet_data(name=task_name,
                            options=options,
                            main_table_name=name,  # 设置交互分析结果目录名
                            module_type="workflow",
                            to_file=to_files,
                            project_sn=project_sn,
                            task_id=data.task_id)

        # 运行workflow 并传回参数
        task_info = super(OnlyTestAction, self).POST()
        task_info['content'] = {
            'ids': {
                'id': str(main_id),
                'name': name
                }
        }
        if 'group_id' in data and str(data.group_id).lower() != 'all':
            _ = self.medical_transcriptome.update_group_is_use(data.task_id, data.group_id)
        if 'control_id' in data:
            _ = self.medical_transcriptome.update_group_compare_is_use(data.task_id, data.control_id)
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
        cmd += "s/medical_transcriptome/diffexp "
        cmd += "-b http://192.168.12.101:9090 "
        args = dict(
            task_id="tsg_33912",
            task_type="2",
            submit_location="diffexp",
            exp_id='5cb9e71b17b2bf450db8efa9',
            group_id='5cb9e50917b2bf450da73f28',
            level='G',
            group_dict=r'{"A1":["A1_1", "A1_2"],"A2": [ "A2_1", "A2_2"], "B1": [ "B1_1", "B1_2"], "B2": [ "B2_1", "B2_2"]}'.replace('"', '\\"'),
            type='ref',
            control_id="5cb9e50917b2bf450da73f29",
            diff_method='DESeq2',
            stat_type='padjust',
            stat_cutoff='0.05',
            correct_method="BY",
            fc='1.5',
            tpm_filter_threshold='0',
        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(str(x) for x in arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTest(TestFunction('test_this'))
    unittest.TextTestRunner(verbosity=2).run(suite)

# -*- coding: utf-8 -*-
import web
import json
import datetime
from collections import OrderedDict
from mainapp.controllers.project.lnc_rna_controller import LncRnaController
from mbio.api.to_file.lnc_rna import *
from mainapp.libs.signature import check_sig
import os
import unittest
from bson.objectid import ObjectId


class DiffexpAction(LncRnaController):
    def __init__(self):
        super(DiffexpAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        basic_args = ["task_id", "submit_location", 'task_type']
        basic_args += ['group_id', 'group_dict', 'control_id',
                       'exp_id', 'exp_level', 'diff_method',
                       'stat_type', 'stat_cutoff', 'fc', 'type','rna_type']#新增rna_type
        # check arg
        for arg in basic_args:
            if not hasattr(data, arg):
                variables = []
                variables.append(arg)
                info = {'success': False, 'info': "Lack argument:   %s" % arg, 'code': 'C2900301', 'variables': variables}
                return json.dumps(info)
        exp_info = self.lnc_rna.get_exp_params_info(data.exp_id, data.task_id)
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
                success.append('If some groups have biological replicates and some groups no, please select DEGseq or edgeR')
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
            fc=data.fc,
            stat_type=data.stat_type,
            stat_cutoff=data.stat_cutoff,
            # quant_method=quant_method,
            diff_method=data.diff_method,
            type=data.type,
            rna_type=data.rna_type
        )
        # for special args
        if hasattr(data, 'correct_method'):
            correct_method = data.correct_method
            params.update(dict(correct_method=correct_method))
        else:
            correct_method = "BH"

        if hasattr(data, 'filter_method'):
            params.update(dict(filter_method=data.filter_method))
            params.update(dict(tpm_filter_threshold=data.tpm_filter_threshold))

        params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        name = "Diff" + '_' + exp_level + '_' + quant_method + '_' +data.rna_type+"_"+data.diff_method + '_'
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
            exp_id=ObjectId(data.exp_id),
            params=params,
            rna_type=data.rna_type,
            status="start",
        )
        main_id = self.lnc_rna.insert_main_table('sg_diff', main_info)

        # prepare option for workflow
        options = {
            "exp_matrix": data.exp_id,
            "group_dict": json.dumps(group_dict),
            "diff_main_id": str(main_id),
            "method": data.diff_method,
            "count": self.lnc_rna.get_count_path(data.exp_id),
            "seq_type_path":self.lnc_rna.get_seq_type_path(data.task_id,exp_level),
            "fc": float(data.fc),
            "group": json.dumps(group_dict),
            "cmp": data.control_id,
            "pvalue": float(data.stat_cutoff),
            "rna_type": data.rna_type,
            "pvalue_padjust": data.stat_type,
            "padjust_way": correct_method,
            "update_info": json.dumps({str(main_id): "sg_diff"}),
            "type": data.type,
            "task_id": data.task_id,
            "exp_level": exp_level,
            "exp_type": exp_type,
        }

        if hasattr(data, 'filter_method') and data.filter_method != "None":
            options.update(dict(filter_method=data.filter_method))
            options.update(dict(tpm_filter_threshold=data.tpm_filter_threshold))

        # prepare to file
        to_files = ["lnc_rna.export_exp_matrix(exp_matrix)",
                    "lnc_rna.export_group(group)",
                    "lnc_rna.export_compare(cmp)",
                    ]

        # 把参数交给workflow运行相应的tool， 其中to_file用于准备tool的输入文件
        task_name = 'lnc_rna.report.diffexp'
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
            _ = self.lnc_rna.update_group_is_use(data.task_id, data.group_id)
        if 'control_id' in data:
            _ = self.lnc_rna.update_group_compare_is_use(data.task_id, data.control_id)
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
        cmd += "s/lnc_rna/diffexp "
        cmd += "-b http://192.168.12.101:9090 "
        args = dict(
            task_id="tsg_33915",
            task_type="2",
            submit_location="diffexp",
            exp_id='5cbd57f417b2bf5efcc410e5',
            group_id='5cbd4ccb17b2bf5efca74408',
            exp_level='trans',
            group_dict=r'{"Con":["Con1","Con2","Con3","Con4","Con5"],"Vit": ["Vit1","Vit2","Vit3","Vit4","Vit5"]}'.replace('"', '\\"'),
            type='all',
            control_id="5cbd4ccb17b2bf5efca74409",
            diff_method='edgeR',
            stat_type='padjust',
            stat_cutoff='0.05',
            fc='1.5',
            rna_type="mRNA+lncRNA"
        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(str(x) for x in arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTest(TestFunction('test_this'))
    unittest.TextTestRunner(verbosity=2).run(suite)

# -*- coding: utf-8 -*-
import web
import json
import datetime
from collections import OrderedDict
from mainapp.controllers.project.itraq_and_tmt_controller import ItraqTmtController
from mbio.api.to_file.itraq_tmt import *
from mainapp.libs.signature import check_sig
import unittest
import os
from bson.objectid import ObjectId


class DiffexpAction(ItraqTmtController):
    def __init__(self):
        super(DiffexpAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        basic_args = ["task_id", "submit_location", 'task_type']
        basic_args += ['group_id', 'group_dict', 'control_id', 'fc_up',
        'fc_down', 'pvalue', 'diff_method', 'type', 'sig_type']
        # type字段确认"origin"，"latest"
        # check arg
        for arg in basic_args:
            if not hasattr(data, arg):
                info = {'success': False, 'info': "Lack argument: %s", "variables":[arg], "code" : "C1900301"}
                return json.dumps(info)
            if arg.lower() == "null":
                info = {'success': False, 'info': "%s : is null or NULL", "variables":[arg], "code" : "C1900302"}
                return json.dumps(info)
        exp_info = self.itraq_tmt.get_exp_params_info_new(task_id=data.task_id, type="ratio")
        project_sn = exp_info["project_sn"]
        task_id = data.task_id
        group_dict = json.loads(data.group_dict, object_pairs_hook=OrderedDict)
        if str(data.group_id).lower() == 'all':
            samples = group_dict['all']
            group_dict = OrderedDict([(x, [x]) for x in samples])
        specimen_group_compare = self.itraq_tmt.get_control_group(task_id=data.task_id,control_id=ObjectId(data.control_id))["compare_names"]
        group_compare = json.loads(specimen_group_compare)
        group_compare_split = [x.split("|") for x in group_compare]
        compare_names = set([item for sublist in group_compare_split for item in sublist])
        print(compare_names)
        print(99999999)
        len_list = list()
        for i in compare_names:
            len_sample = len(group_dict[i])
            len_list.append(len_sample)
        success = list()
        if len_list.count(1) >= 2:
            if data.diff_method == "student" or data.diff_method == "welch":
                success.append("student's t-test and Welch's T Test are not suitable to deal with single sample and single sample comparison, please select Chisq and Fisher's Exact instead.")
        if len_list.count(1) >= 1:
            if data.diff_method == "welch":
                success.append("Welch's T Test is not suitable to deal with single sample and multiple samples comparison, please select student's t-test instead.")
        if len_list.count(1) != len(compare_names):
            if data.diff_method == "chi" or data.diff_method == "fisher":
                success.append("Deal with multiple samples and multiple samples comparison，we recommend student's t-test and Welch's T Test")
        if success:
            info = {'success': False, 'info': success[0], "code" : "C1900303"}
            return json.dumps(info)

        # create main table record
        params = dict(
            task_id=data.task_id,
            submit_location=data.submit_location,
            task_type=int(data.task_type),
            group_id=data.group_id,
            control_id=data.control_id,
            group_dict=json.loads(data.group_dict, object_pairs_hook=OrderedDict),
            fc_up=data.fc_up,
            fc_down=data.fc_down,
            pvalue=data.pvalue,
            diff_method=data.diff_method,
            type=data.type,
        )
        # 卡方检验的页面没有单双尾检验
        if hasattr(data, 'correct_method'):
            correct_method = data.correct_method
            params.update(dict(correct_method=correct_method))
        if hasattr(data, 'padjust_way'):
            params.update(dict(padjust_way=data.padjust_way))
        if hasattr(data, 'sig_type'):
            params.update(dict(sig_type=data.sig_type))
        params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        name = "Diff" + '_'  + data.diff_method + '_'
        time_now = datetime.datetime.now()
        name += time_now.strftime("%Y%m%d_%H%M%S")
        if type(params) == dict:
            params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        result_dir = self.itraq_tmt.get_anno_info(data.type, data.task_id)["result_dir"]
        if not result_dir.endswith('/'):
            result_dir += '/'
        protein_sliced = self.itraq_tmt.get_task_info(data.task_id)["protein_sliced"]
        main_info = dict(
            project_sn=project_sn,
            task_id=task_id,
            name=name,
            created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
            desc='differential analysis main table',
            group_dict=group_dict,
            group_id=data.group_id,
            params=params,
            status="start",
        )
        #　self.itraq_tmt来自mainapp/controllers/project
        # /itraq_and_tmt_controller.py的from
        # mainapp.models.mongo.itraq_and_tmt import ItraqTmt然后init里面写了self.itraq_tmt = ItraqTmt(bind_object=bind_object)
        # 这里就是接口直接调用函数，而to_file里面的函数是接口造字符串，传递给workflow才使用函数调用
        main_id = self.itraq_tmt.insert_main_table('sg_diff', main_info)

        # prepare option for workflow
        options = {
            "exp_matrix": data.task_id,
            "group_dict": json.dumps(group_dict),
            "fc_up": float(data.fc_up),
            "fc_down": float(data.fc_down),
            "method_type": data.diff_method,
            "group": json.dumps(group_dict),
            "cmp": data.control_id,
            "control_id": data.control_id,
            "result_dir":result_dir,
            "protein_sliced":protein_sliced,
            "diff_main_id": str(main_id),
            "sig_type": data.sig_type,
            "pvalue": data.pvalue,
            "update_info": json.dumps({str(main_id): "sg_diff"}), # to update sg_status
        }
        if hasattr(data, 'correct_method'):
            correct_method = data.correct_method
            options.update(dict(correct_method=correct_method))
        if hasattr(data, 'padjust_way'):
            options.update(dict(padjust_way=data.padjust_way))
        # prepare to file
        to_files = ["itraq_tmt.export_exp_matrix(exp_matrix)",
                    "itraq_tmt.export_group(group)",
                    "itraq_tmt.export_compare(cmp)"]

        # 把参数交给workflow运行相应的tool， 其中to_file用于准备tool的输入文件
        task_name = 'itraq_and_tmt.report.diffexp'
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
        cmd += "s/itraq_and_tmt/diffexp "
        cmd += "-b http://bcl.tsg.com "
        args = dict(
            task_id="tsg_249300",
            task_type="2",
            submit_location="diffexp",
            control_id="5ff639cf17b2bf6c30638f76",
            correct_method="two.sided",
            diff_method="student",
            fc_down="0.83",
            fc_up="1.2",
            group_dict=json.dumps({"C": ["C1", "C2", "C3"], "E": ["E1", "E2", "E3"]}).replace('"', '\\"'),
            group_id="5ff639cf17b2bf6c30638f75",
            sig_type='padjust',
            pvalue="0.05",
            padjust_way='3',
            type="origin"

        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)


if __name__ == '__main__':
    unittest.main()

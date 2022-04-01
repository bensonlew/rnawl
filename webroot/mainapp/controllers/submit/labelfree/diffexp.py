# -*- coding: utf-8 -*-
import web
import json
import datetime
from collections import OrderedDict
from mainapp.controllers.project.labelfree_controller import LabelfreeController
from mbio.api.to_file.labelfree import *
from mainapp.libs.signature import check_sig
import unittest
import os


class DiffexpAction(LabelfreeController):
    def __init__(self):
        super(DiffexpAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        basic_args = ["task_id", "submit_location", 'task_type']
        basic_args += ['group_id', 'group_dict', 'control_id', 'fc_up',
                       'fc_down', 'pvalue', 'diff_method', 'type']
        # type字段确认"origin"，"latest"
        # check arg
        for arg in basic_args:
            if not hasattr(data, arg):
                info = {'success': False, 'info': "Lack argument: {}".format(arg)}
                return json.dumps(info)
            if arg.lower() == "null":
                info = {'success': False, 'info': "{} : is null or NULL".format(arg)}
                return json.dumps(info)
        exp_info = self.labelfree.get_exp_params_info_new(task_id=data.task_id, type="ratio")
        project_sn = exp_info["project_sn"]
        task_id = data.task_id
        group_dict = json.loads(data.group_dict, object_pairs_hook=OrderedDict)
        if str(data.group_id).lower() == 'all':
            samples = group_dict['all']
            group_dict = OrderedDict([(x, [x]) for x in samples])
        specimen_group_compare = \
        self.labelfree.get_control_group(task_id=data.task_id, control_id=ObjectId(data.control_id))["compare_names"]
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
                success.append(
                    "student's t-test" + "和" "Welch's T Test" + "不适合处理单样本和单样本比较, 请选择Chisq或Fisher's Exact,并且确认你的每组只有一个样本")
        if len_list.count(1) >= 1:
            if data.diff_method == "welch":
                success.append("Welch's T Test" + "不适合处理组内同时有单样本和多样本的, 此时请选择student's t-test")
        if len_list.count(1) != len(compare_names):
            if data.diff_method == "chi" or data.diff_method == "fisher":
                success.append("只涉及包含多样本的组与组间比较时，我们只推荐 " + "student's ""t-test" + "和" + "Welch's T Test")
        if success:
            info = {'success': False, 'info': success[0]}
            return json.dumps(info)

        # create main table record
        params = dict(
            task_id=data.task_id,
            submit_location=data.submit_location,
            task_type=int(data.task_type),
            group_id=data.group_id,
            group_dict=json.loads(data.group_dict, object_pairs_hook=OrderedDict),
            fc_up=data.fc_up,
            fc_down=data.fc_down,
            pvalue=data.pvalue,
            # correct_method=data.correct_method,
            diff_method=data.diff_method,
            control_id=data.control_id,
            type=data.type,
        )
        if hasattr(data, 'cutoffs'):
            params.update({"cutoffs":data.cutoffs})
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
        result_dir = self.labelfree.get_anno_info(data.type, data.task_id)["result_dir"]
        if not result_dir.endswith('/'):
            result_dir += '/'
        protein_sliced = self.labelfree.get_task_info(data.task_id)["protein_sliced"]
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
        #　self.labelfree来自mainapp/controllers/project
        # /labelfree_controller.py的from
        # mainapp.models.mongo.labelfree import ItraqTmt然后init里面写了self.labelfree = ItraqTmt(bind_object=bind_object)
        # 这里就是接口直接调用函数，而to_file里面的函数是接口造字符串，传递给workflow才使用函数调用
        main_id = self.labelfree.insert_main_table('sg_diff', main_info)


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
            # "sig_type": data.sig_type,
            "pvalue": data.pvalue,
            "update_info": json.dumps({str(main_id): "sg_diff"})  # to update sg_status
        }
        if hasattr(data, 'sig_type'):
            options.update({"sig_type":data.sig_type})
        if hasattr(data, 'cutoffs'):
            options.update({"cutoffs":data.cutoffs})
        if hasattr(data, 'correct_method'):
            correct_method = data.correct_method
            options.update(dict(correct_method=correct_method))
        if hasattr(data, 'padjust_way'):
            options.update(dict(padjust_way=data.padjust_way))
        # prepare to file
        to_files = ["labelfree.export_exp_matrix(exp_matrix)",
                    "labelfree.export_group(group)",
                    "labelfree.export_compare(cmp)"]

        # 把参数交给workflow运行相应的tool， 其中to_file用于准备tool的输入文件
        task_name = 'labelfree.report.diffexp'
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
        cmd += "s/labelfree/diffexp "
        cmd += "-b http://bcl.tsg.com "
        args = dict(
            task_id="tsg_34739",
            task_type="2",
            submit_location="diffexp",
            group_id="5d22a82717b2bf1c689db738",
            group_dict=r'{"F":["F_1","F_2","F_3"],"S":["S_1","S_2","S_3"]}'.replace(
                '"', '\\"'),
            control_id="5d22a82717b2bf1c689db739",
            fc_up="1.2",
            fc_down="0.83",
            pvalue="0.05",
            diff_method="student",
            correct_method="two.sided",
            type="origin"

        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)


if __name__ == '__main__':
    unittest.main()

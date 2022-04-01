# -*- coding: utf-8 -*-
import web
import json
import datetime
from collections import OrderedDict
from mainapp.controllers.project.ref_rna_v2_controller import RefRnaV2Controller
from mbio.api.to_file.ref_rna_v2 import *
from mainapp.libs.signature import check_sig
import unittest
import os


class SangerChartAction(RefRnaV2Controller):
    def __init__(self):
        super(SangerChartAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        print(data)
        basic_args = ["task_id", "submit_location", 'task_type']
        basic_args += ['url','content']
        # check arg
        for arg in basic_args:
            if not hasattr(data, arg):
                variables = []
                variables.append(arg)
                info = {'success': False, 'info': "Lack argument: %s" % arg, 'code': 'C2900701', 'variables': variables}
                return json.dumps(info)
        exp_info = self.ref_rna_v2.get_task_info(task_id=data.task_id)
        project_sn = exp_info["project_sn"]
        task_id = data.task_id
        if not hasattr(data, 'user_name'):
            user_name = 'sanger_chart@majorbio.com'
            password = 'sanger1234'
            userpass = user_name + '____username_password____' + password
        else:
            user_name = data.user_name
            password = data.password
            userpass = user_name + '____username_password____' + password
        # create main table record
        params = dict(
            task_id=task_id,
            submit_location=data.submit_location,
            task_type=int(data.task_type),
            url=data.url,
            user=user_name,
            password=password,
            content=data.content,
        )
        params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        name = "SangerChart" + '_'
        time_now = datetime.datetime.now()
        name += time_now.strftime("%Y%m%d_%H%M%S")
        if type(params) == dict:
            params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        main_info = dict(
            project_sn=project_sn,
            task_id=task_id,
            version="v1",
            name=name,
            created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
            desc='Sanger Chart main table',
            params=params,
            status="start"
        )
        main_id = self.ref_rna_v2.insert_main_table('sanger_chart', main_info)
        new_task_id = self.ref_rna_v2.get_new_id(data.task_id)
        main_table_data = {'run_id': new_task_id}

        # prepare option for workflow
        options = {
            "url":data.url,
            "userpass":userpass,
            "mode":data.content,
            "chart_main_id": str(main_id),
            "update_info": json.dumps({str(main_id): "sanger_chart"})  # to update sg_status
        }

        # 把参数交给workflow运行相应的tool， 其中to_file用于准备tool的输入文件
        task_name = 'ref_rna_v2.report.sanger_chart'
        self.set_sheet_data(name=task_name,
                            options=options,
                            main_table_name=name,  # 设置交互分析结果目录名
                            module_type="workflow",
                            project_sn=project_sn,
                            new_task_id=new_task_id,
                            task_id=data.task_id)

        # 运行workflow 并传回参数
        task_info = super(SangerChartAction, self).POST()
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
        cmd = 'python /mnt/ilustre/users/isanger/sanger_bioinfo/bin/webapitest_new.py '
        cmd += 'post '
        cmd += "-fr no "
        cmd += '-c {} '.format("client03")
        cmd += "s/ref_rna_v2/sanger_chart "
        cmd += "-b http://bcl.i-sanger.com "
        cmd += "-dbversion 1 "
        args = dict(
            task_id="majorbio_314991",
            task_type="2",
            submit_location="sangerchart",
            url='http://report.sanger.com/refrna/genesetvenn/task_id/majorbio_314991.html',
            user_name='sgtest@majorbio.com',
            password='test01',
            content='single'
        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(str(x) for x in arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTest(TestFunction('test_this'))
    unittest.TextTestRunner(verbosity=2).run(suite)

# -*- coding: utf-8 -*-
import web
import json
import datetime
from collections import OrderedDict
from mainapp.controllers.project.dia_controller import DiaController
from mbio.api.to_file.dia import *
from mainapp.libs.signature import check_sig
import unittest
import os

class ExpPcaAction(DiaController):
    def __init__(self):
        super(ExpPcaAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        basic_args = ["task_id", "submit_location", 'task_type']
        basic_args += ['group_id', 'group_dict']
        basic_args += ['express_id']
        # check arg
        for arg in basic_args:
            if not hasattr(data, arg):
                info = {'success': False, 'info': "Lack argument: {}".format(arg)}
                return json.dumps(info)
            if not hasattr(data, arg):
                info = {'success': False, 'info': "{} : is null or NULL".format(arg)}
                return json.dumps(info)
        # exp_info = self.dia.get_exp_params_info_new(task_id=data.task_id, type="ratio")
        exp_info = self.dia.get_main_info(main_id=data.express_id, collection_name="sg_express", task_id=data.task_id)
        project_sn = exp_info["project_sn"]
        task_id = data.task_id
        group_dict = json.loads(data.group_dict, object_pairs_hook=OrderedDict)
        # create main table record
        params = dict(
            task_id=task_id,
            submit_location=data.submit_location,
            task_type=int(data.task_type),
            express_id=data.express_id,
            group_id=data.group_id,
            group_dict=group_dict,
        )
        params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        name = "SamPCA" + '_'
        time_now = datetime.datetime.now()
        name += time_now.strftime("%Y%m%d_%H%M%S")
        if type(params) == dict:
            params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        main_info = dict(
            project_sn=project_sn,
            task_id=task_id,
            name=name,
            created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
            desc='PCA main table',
            params=params,
            status="start",
            version='v3',
        )
        main_id = self.dia.insert_main_table('sg_express_pca', main_info)

        # prepare option for workflow
        if str(data.group_id).lower() == 'all':
            samples = group_dict['all']
            group_dict = OrderedDict([(x, [x]) for x in samples])
        options = {
            "exp_matrix": data.express_id,
            "group_dict": json.dumps(group_dict),
            "group": json.dumps(group_dict),
            "pca_main_id": str(main_id),
            "update_info": json.dumps({str(main_id): "sg_express_pca"})  # to update sg_status
        }
        # prepare to file
        to_files = ["dia.export_exp_matrix_new(exp_matrix)",
                    "dia.export_group(group)"]

        # 把参数交给workflow运行相应的tool， 其中to_file用于准备tool的输入文件
        task_name = 'dia_v3.report.exp_pca'
        self.set_sheet_data(name=task_name,
                            options=options,
                            main_table_name=name,  # 设置交互分析结果目录名
                            module_type="workflow",
                            to_file=to_files,
                            project_sn=project_sn,
                            task_id=data.task_id)

        # 运行workflow 并传回参数
        task_info = super(ExpPcaAction, self).POST()
        task_info['content'] = {
            'ids': {
                'id': str(main_id),
                'name': name
                }
        }
        # task_info['group_dict'] = group_dict

        ##
        if 'group_id' in data and str(data.group_id).lower() != 'all':
            _ = self.dia.update_group_is_use(data.task_id, data.group_id)
        if 'control_id' in data:
            _ = self.dia.update_group_compare_is_use(data.task_id, data.control_id)
        ##
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
        cmd += "s/dia_v3/exp_pca "
        cmd += '-b http://bcl.tsg.com '
        args = dict(
            task_id="dia_test",
            task_type="2",
            submit_location="expresspca",
            express_id='5f7763e317b2bf57d7ed373b',
            group_id="5f7763e317b2bf57d7ed3742",
            group_dict=json.dumps({"C_3": ["C_3_1", "C_3_2", "C_3_3", "C_3_4", "C_3_5"],
                                   "R_3": ["R_3_1", "R_3_2", "R_3_3", "R_3_4", "R_3_5"]}).replace('"', '\\"'),
        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)


if __name__ == '__main__':
    unittest.main()

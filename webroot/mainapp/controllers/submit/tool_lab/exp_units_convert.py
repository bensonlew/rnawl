# -*- coding: utf-8 -*-
# __author__ = 'shicaiping'
import codecs
import datetime
import json
import os
import unittest

import web

from mainapp.controllers.project.tool_lab_controller import ToolLabController
from mainapp.libs.signature import check_sig


class ExpUnitsConvertAction(ToolLabController):
    def __init__(self):
        super(ExpUnitsConvertAction, self).__init__(instant=False)
        self.basic_args = ["task_id", "exp_matrix", "from_unit", "to_unit", "float_num", "gene_length", "intersect", "task_type", "submit_location"]
        self.input_data = web.input()

    @check_sig
    def POST(self):
        # check params
        try:
            check = self.check_params()
            if check is not True:
                return check
        except Exception, e:
            info = {'success': False, 'info': '{}'.format(e)}
            return json.dumps(info)

        # create main table
        try:
            params = self.pack_params()
            name = "ExpUnitsConvert" + '_' + self.input_data.from_unit.lower() + "2" + self.input_data.to_unit.lower() + '_'
            time_now = datetime.datetime.now()
            name += time_now.strftime("%Y%m%d_%H%M%S")
            main_info = dict(
                project_sn=self.input_data.project_sn,
                task_id=self.input_data.task_id,
                name=name,
                version="v1",
                created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
                desc=self.input_data.from_unit.lower() + "2" + self.input_data.to_unit.lower(),
                params=params,
                status="start"
            )
            main_id = self.tool_lab.insert_main_table('exp_units_convert', main_info)
        except Exception, e:
            info = {'success': False, 'info': '{}'.format(e)}
            return json.dumps(info)

        # prepare option for workflow
        try:
            options = {
                "exp_matrix": self.input_data.exp_matrix,
                "convert_type": self.input_data.from_unit.lower() + "2" + self.input_data.to_unit.lower(),
                "gene_length": self.input_data.gene_length,
                "float_num": self.input_data.float_num,
                "intersect": self.input_data.intersect,
                "update_info": json.dumps({str(main_id): "exp_units_convert"})  # to update sg_status
            }

            # 把参数交给workflow运行相应的tool
            task_name = 'tool_lab.exp_units_convert'
            self.set_sheet_data(name=task_name,
                                options=options,
                                main_table_name=name,  # 设置交互分析结果目录名
                                module_type="workflow",
                                project_sn=self.input_data.project_sn,
                                task_id=self.input_data.task_id)
        except Exception, e:
            info = {'success': False, 'info': '{}'.format(e)}
            return json.dumps(info)

        # 运行workflow 并传回参数
        task_info = super(ExpUnitsConvertAction, self).POST()
        task_info['content'] = {
            'ids': {
                'id': str(main_id),
                'name': name
            }
        }
        return json.dumps(task_info)

    def check_params(self):
        for arg in self.basic_args:
            if not hasattr(self.input_data, arg):
                variables = []
                variables.append(arg)
                info = {'success': False, 'info': "Lack argument: %s" % arg, 'code': 'C2901801', 'variables': variables}
                return json.dumps(info)
        convert = self.input_data.from_unit.lower() + "2" + self.input_data.to_unit.lower()
        if convert not in ['count2cpm', 'fpkm2tpm', 'count2tpm', 'count2fpkm', 'cpm2fpkm', 'cpm2tpm', 'fpkm2cpm']:
            info = {'success': False, 'info': "%s to %s not supported!" % (self.input_data.from_unit, self.input_data.to_unit)}
            return json.dumps(info)
        if convert not in ['count2cpm', 'fpkm2tpm']:
            if not self.input_data.gene_length or self.input_data.gene_length == "null":
                info = {'success': False,
                        'info': "%s to %s should provie gene length file!" % (self.input_data.from_unit, self.input_data.to_unit)}
                return json.dumps(info)
        return True

    def pack_params(self):
        input_data_dict = dict(self.input_data)
        params_dict = dict()
        args = self.basic_args
        for each in args:
            if each == "task_type":
                params_dict[each] = int(input_data_dict[each])
            else:
                params_dict[each] = input_data_dict[each]
        return json.dumps(params_dict, sort_keys=True, separators=(',', ':'))

class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """

    def test_this(self):
        cmd = 'python /mnt/ilustre/users/sanger-dev/biocluster/bin/webapitest.py '
        cmd += 'post '
        cmd += "-fr no "
        cmd += '-c {} '.format("client03")
        cmd += "s/tool_lab/exp_units_convert "
        cmd += "-b http://bcl.tsg.com "
        args = dict(
            task_id="exp_units_convert",
            project_sn="exp_units_convert",
            task_type="2",
            submit_location="exp_units_convert",
            name="exp_units_convert",
            exp_matrix="/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/de_tools/known_seqs_count.matrix",
            gene_length="/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/de_tools/gene_length",
            from_unit='count',
            to_unit='tpm',
            float_num='4',
            intersect='True',
        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)


if __name__ == '__main__':
    unittest.main()

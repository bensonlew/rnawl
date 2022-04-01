# -*- coding: utf-8 -*-
# __author__ = 'fwy'

import codecs
import datetime
import json
import os
import unittest

import web

from mainapp.controllers.project.tool_lab_controller import ToolLabController
from mainapp.libs.signature import check_sig


class SeqFilterAction(ToolLabController):
    def __init__(self):
        super(SeqFilterAction, self).__init__(instant=False)
        self.basic_args = ["task_id", "is_upload", "min_len", "max_len","fasta_file","fasta_seq"]
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
        params = self.pack_params()
        name = "SeqFilter" + '_' + "min"+self.input_data.min_len + "_" + "max"+self.input_data.max_len + '_'
        time_now = datetime.datetime.now()
        name += time_now.strftime("%Y%m%d_%H%M%S")
        main_info = dict(
            project_sn=self.input_data.project_sn,
            task_id=self.input_data.task_id,
            name=name,
            version="v1",
            created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
            desc="SeqFilter",
            params=params,
            status="start"
        )
        main_id = self.tool_lab.insert_main_table('seq_filter', main_info)


        # prepare option for workflow
        # try:
        options = {
            "is_upload": self.input_data.is_upload,
            "max_len": self.input_data.max_len,
            "min_len": self.input_data.min_len,
            "fasta_seq": self.input_data.fasta_seq,
            "update_info": json.dumps({str(main_id): "seq_filter",})  # to update sg_status
        }
        if self.input_data.is_upload == "yes":
            options.update({"fasta_file": self.input_data.fasta_file})

        # 把参数交给workflow运行相应的tool
        task_name = 'tool_lab.seq_filter'
        self.set_sheet_data(name=task_name,
                            options=options,
                            main_table_name=name,  # 设置交互分析结果目录名
                            module_type="workflow",
                            project_sn=self.input_data.project_sn,
                            task_id=self.input_data.task_id)
        # except Exception, e:
        #     info = {'success': False, 'info': '{}'.format(e)}
        #     return json.dumps(info)

        # 运行workflow 并传回参数
        task_info = super(SeqFilterAction, self).POST()
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
        try:
            min_len=int(self.input_data.min_len)
        except:
            info = {'success': False,
                    'info': "最小长度必须输入整数"}
            return json.dumps(info)
        try:
            max_len = int(self.input_data.min_len)
        except:
            info = {'success': False,
                    'info': "最大长度必须输入整数"}
            return json.dumps(info)
        if min_len >= max_len:
            info = {'success': False,
                    'info': "最小长度必须小于最大长度"}
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
        cmd += "s/tool_lab/seq_filter "
        cmd += "-b http://bcl.tsg.com "
        args = dict(
            task_id="seq_filter",
            project_sn="seq_filter",
            task_type="2",
            submit_location="seq_filter",
            # name="seq_filter",
            fasta_file="/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/tool_lab/seq_filter/input/example.fasta",
            min_len='500.1',
            max_len='10000',
            fasta_seq="",
            is_upload="yes"
        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)


if __name__ == '__main__':
    unittest.main()

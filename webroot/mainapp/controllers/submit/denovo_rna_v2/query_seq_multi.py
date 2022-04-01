# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'

import os
import web
import json
import datetime
from mainapp.controllers.project.denovo_rna_v2_controller import DenovoRnaV2Controller
from mainapp.libs.signature import check_sig
import sqlite3
import unittest
import re

class QuerySeqMultiAction(DenovoRnaV2Controller):
    def __init__(self):
        super(QuerySeqMultiAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        # check
        for arg in ['seq_id', 'seq_type', 'task_id']:
            if not hasattr(data, arg):
                var = []
                var.append(arg)
                info = {'success': False, 'info': "Lack argument: %s" , "variables":[arg], "code" : "C1603001"}
                return json.dumps(info)

        if data['seq_type'] not in ['T', 'G', 'P']:
            info = {'success': False, 'info': 'unexpected seq_type', "code" : "C1603002"}
            return json.dumps(info)

        task_info = self.denovo_rna_v2.get_task_info(data.task_id)

        params = dict(
            task_id=data.task_id,
            task_type=int(data.task_type),
            submit_location=data.submit_location,
            seq_id=data.seq_id,
            seq_type=data.seq_type
        )
        print "params is {}".format(params)
        if type(params) == dict:
            params = json.dumps(params, sort_keys=True, separators=(',', ':'))

        name = "Download_sequence_"
        time_now = datetime.datetime.now()
        name += time_now.strftime("%Y%m%d_%H%M%S")
        main_table_name = "download_data"
        main_info = dict(
            project_sn=task_info["project_sn"],
            task_id=task_info["task_id"],
            name=name,
            created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
            desc='download sequence',
            params=params,
            status="start",
        )
        main_id = self.denovo_rna_v2.insert_main_table('sg_query_seq_multi', main_info)

        options = {
            "seq_id": data.seq_id,
            "seq_type": data.seq_type,
            "seq_db": task_info["seq_db"],
            "task_id": task_info["task_id"],
            "main_id": str(main_id)

        }

        # 不跟新运行记录
        # "update_info": json.dumps({str(main_id): "sg_query_seq_multi"})

        to_files = []
        project_sn = task_info["project_sn"]
        task_name = 'denovo_rna_v2.report.query_seq_multi'
        self.set_sheet_data(name=task_name,
                            options=options,
                            main_table_name=name,  # 设置交互分析结果目录名
                            module_type="workflow",
                            to_file=to_files,
                            project_sn=project_sn,
                            task_id=data.task_id)

        # 运行workflow 并传回参数
        task_info = super(QuerySeqMultiAction, self).POST()
        task_info['content'] = {
            'ids': {
                'id': str(main_id),
                'name': name
                }
        }
        if 'group_id' in data and str(data.group_id).lower() != 'all':
            _ = self.denovo_rna_v2.update_group_is_use(data.task_id, data.group_id)
        if 'control_id' in data:
            _ = self.denovo_rna_v2.update_group_compare_is_use(data.task_id, data.control_id)
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
        cmd += "s/denovo_rna_v2/query_seq_multi "
        cmd += "-b http://192.168.12.101:51232 "
        args = dict(
            task_id="denovo_rna_v2_upgrade",
            task_type="2",
            submit_location="query_seq_multi",
            seq_type="T",
            seq_id="TRINITY_DN10012_c0_g1_i1,TRINITY_DN10051_c0_g1_i1,TRINITY_DN10062_c0_g1_i1,TRINITY_DN10178_c0_g1_i1,TRINITY_DN1021_c0_g1_i1,TRINITY_DN10311_c0_g1_i1"

        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(str(x) for x in arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTest(TestFunction('test_this'))
    unittest.TextTestRunner(verbosity=2).run(suite)

# -*- coding: utf-8 -*-
# __author__ = 'shicaiping'

import codecs
import datetime
import os
import unittest

import web
import json
import re

from mainapp.controllers.project.prok_rna_controller import ProkRNAController
from mainapp.libs.signature import check_sig
from mbio.api.to_file.prok_rna import *


class GenesetBatchAction(ProkRNAController):
    def __init__(self):
        super(GenesetBatchAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        input_data = web.input()
        print(input_data)
        basic_args = ["task_id", "submit_location", "diff_id", "task_type"]
        # check arg
        for arg in basic_args:
            if not hasattr(input_data, arg):
                variables = list()
                variables.append(arg)
                info = {'success': False, 'info': "Lack argument: %s" % arg, 'code': 'C2901801', 'variables': variables}
                return json.dumps(info)
            if arg.lower() == "null":
                variables = list()
                variables.append(arg)
                info = {'success': False, 'info': "%s : is null or NULL" % arg, 'code': 'C2901802',
                        'variables': variables}
                return json.dumps(info)

        diff_info = self.prok_rna.get_main_info(input_data.diff_id, "sg_diff", input_data.task_id)
        diff_name = diff_info["name"]
        project_sn = diff_info["project_sn"]
        geneset_type = diff_info["exp_level"]
        cmp_detail = diff_info["cmp_detail"]
        if diff_info["status"] != "end":
            info = {'success': False, 'info': "运行没有结束或已经失败，无法批量创建基因集"}
            return json.dumps(info)
        params = json.loads(diff_info["params"])
        flag = int(self.prok_rna.get_geneset_batch_flag(input_data.task_id))
        if "is_batch" in params and params["is_batch"] == "True":
            geneset_suffix = "_batch_" + str(flag)
        else:
            geneset_suffix = "_" + str(flag)
        print(geneset_suffix)
        geneset_names = self.prok_rna.get_genesets(input_data.task_id, geneset_type, is_use=1)
        for cmp in cmp_detail:
            if geneset_suffix:
                geneset_name = cmp.replace("|", "_vs_") + "_" + geneset_type + geneset_suffix
            else:
                geneset_name = cmp.replace("|", "_vs_") + "_" + geneset_type
            if geneset_name in geneset_names:
                info = {'success': False, 'info': "基因集{}存在分析记录，请在基因集管理处先删除该基因集后再创建".format(geneset_name)}
                return json.dumps(info)
        # create main table record
        params = dict(
            task_id=input_data.task_id,
            submit_location=input_data.submit_location,
            task_type=int(input_data.task_type),
            diff_id=input_data.diff_id
        )
        params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        time = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        desc = "batch create genesets of table {}".format(diff_name)
        name = "Batch_create_genesets_{}".format(time)
        # prepare main table info
        main_info = dict(
            name=name,
            task_id=input_data.task_id,
            project_sn=project_sn,
            desc=desc,
            status="start",
            params=params,
            flag=flag,
            created_ts=datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
        )
        main_id = self.prok_rna.insert_main_table('sg_geneset_batch', main_info)
        new_task_id = self.prok_rna.get_new_id(input_data.task_id)
        main_table_data = {'run_id': new_task_id}

        # prepare option for workflow
        options = {
            "genesets": input_data.diff_id,
            "geneset_type": geneset_type,
            "geneset_suffix": geneset_suffix,
            "diff_id": input_data.diff_id,
            "task_id": input_data.task_id,
            "main_table_data": main_table_data,
            "update_info": json.dumps({str(main_id): "sg_geneset_batch"})  # to update sg_status
        }
        # prepare to file
        # to_files = ["ref_rna_v2.export_diff_genesets(genesets)"]

        # 把参数交给workflow运行相应的tool， 其中to_file用于准备tool的输入文件
        task_name = 'prok_rna.report.geneset_batch'
        self.set_sheet_data(name=task_name,
                            options=options,
                            main_table_name=name,  # 设置交互分析结果目录名
                            module_type="workflow",
                            # to_file=to_files,
                            project_sn=project_sn,
                            new_task_id=new_task_id,
                            task_id=input_data.task_id)

        # 运行workflow 并传回参数
        task_info = super(GenesetBatchAction, self).POST()
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
        cmd += "s/ref_rna_v3/geneset_batch "
        cmd += "-b http://bcl.tsg.com "
        args = dict(
            task_id="tsg_218497",
            task_type="2",
            submit_location="geneset_batch",
            diff_id='5f4878fc17b2bf305c20ad46'
        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)


if __name__ == '__main__':
    unittest.main()

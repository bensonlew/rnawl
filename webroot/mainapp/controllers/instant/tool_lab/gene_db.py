# -*- coding: utf-8 -*-

import datetime
import json
import os
import unittest
from collections import OrderedDict
import web
from mainapp.controllers.project.tool_lab_controller import ToolLabController
from mainapp.libs.signature import check_sig
import random

class GeneDbAction(ToolLabController):
    def __init__(self):
        super(GeneDbAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        basic_args = ["task_id", "task_type", "submit_location", "task_name", "species", "search_field", "target_field"]
        all_args =  basic_args + ['search_file', "search_string", "file_id"]
        # check arg

        for arg in basic_args:
            if not hasattr(data, arg):
                variables = []
                variables.append(arg)
                info = {'success': False, 'info': "Lack argument: %s" % arg}
                return json.dumps(info)

        task_id = data.task_id
        params = dict()
        for arg in all_args:
            if hasattr(data, arg):
                if arg == "task_type":
                    params.update({arg: int(getattr(data, arg))})
                else:
                    params.update({arg: getattr(data, arg)})

        name = "GeneDb" + '_'
        if data.task_name != "":
            name += data.task_name + "_"
        time_now = datetime.datetime.now()
        name += time_now.strftime("%Y%m%d_%H%M%S")
        if type(params) == dict:
            params = json.dumps(params, sort_keys=True, separators=(',', ':'))

        main_info = dict(
            project_sn=data.project_sn,
            task_id=data.task_id,
            name=name,
            version="v1",
            created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
            desc='gene database id convert',
            params=params,
            status="start"
        )
        main_id = self.tool_lab.insert_main_table('sgdb_id_search', main_info)

        options = dict()
        for arg in all_args:
            if hasattr(data, arg):
                options.update({arg: getattr(data, arg)})

        options.update({
            "main_id": str(main_id)
        })
        options.update({
            "update_info": json.dumps({str(main_id): "sgdb_id_search"})
        })

        to_files = []
        # 把参数交给workflow运行相应的tool， 其中to_file用于准备tool的输入文件
        task_name = 'tool_lab.genedb_search'
        self.set_sheet_data(name=task_name,
                            options=options,
                            main_table_name=name,  # 设置交互分析结果目录名
                            module_type="workflow",
                            to_file=to_files,
                            project_sn=data.project_sn,
                            task_id=data.task_id)

        # 运行workflow 并传回参数
        task_info = super(GeneDbAction, self).POST()
        task_info['content'] = {
            'ids': {
                'id': str(main_id),
                'name': name
            }
        }

        return json.dumps(task_info)


class TestFunction(unittest.TestCase):
    """
    This is test for the workflow. Just run this script to do test.
    """


    def test_this(self):
        cmd = 'python /mnt/ilustre/users/sanger-dev/biocluster/bin/webapitest.py '
        cmd += 'post '
        cmd += "-fr no "
        cmd += '-c {} '.format("client03")
        cmd += "i/tool_lab/gene_db "
        cmd += "-b http://bcl.tsg.com "
        args = dict(
            task_id="test_" + str(random.randint(1, 10000)),
            project_sn="test_123",
            task_type="1",
            task_name="test_123",
            submit_location="gene_db",
            species="hsapiens",
            search_field="entrezgene_id,ensembl_gene_id,uniprot_gn_id,external_gene_name",
            target_field="all",
            search_string="101928721,ENSG00000273252,Q9Y225,RNF24,ENSG00000148634,ENSG00000265590,57019,Q9BYV7,CIAPIN1,SCAF4"

        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)


if __name__ == '__main__':
    unittest.main()

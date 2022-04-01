# -*- coding: utf-8 -*-
# __author__ = 'qindanhua'
import web
import json
import datetime
from mainapp.controllers.project.ref_rna_v2_controller import RefRnaV2Controller
from mbio.api.to_file.ref_rna_v2 import *
from mainapp.libs.signature import check_sig
import os
import unittest

class GenesetKeggAction(RefRnaV2Controller):
    def __init__(self):
        super(GenesetKeggAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        default_argu = ['geneset_id', 'geneset_type', 'submit_location', 'type', 'task_type', 'task_id']
        for argu in default_argu:
            if not hasattr(data, argu):
                variables = []
                variables.append(argu)
                info = {'success': False, 'info': '%s参数缺少!' % argu, 'code': 'C2901601', 'variables': variables}
                return json.dumps(info)

        task_name = 'ref_rna_v2.report.geneset_kegg'
        params_json = {
            "submit_location": data.submit_location,
            "task_type": int(data.task_type),
            "geneset_id": data.geneset_id,
            "geneset_type": data.geneset_type,
            "task_id": data.task_id,
            "type": data.type,
        }
        # 判断传入的基因集id是否存在
        geneset_info = {}
        geneset_documents = []
        if len(data.geneset_id.split(",")) >=3:
            info = {"success": False, "info": "只能选择一个或两个基因集", "code" : "C2901604"}
            return json.dumps(info)
        for geneset in data.geneset_id.split(","):
            geneset_info = self.ref_rna_v2.get_main_info(geneset, 'sg_geneset', data.task_id)
            if not geneset_info:
                info = {"success": False, "info": "geneset不存在，请确认参数是否正确！!", 'code': 'C2901602', 'variables': ''}
                return json.dumps(info)
            geneset_documents.append(geneset_info)
        source = 'non_diff_exp'
        if len(geneset_documents) == 1:
            my_result = geneset_documents[0]
            if 'source' in my_result and my_result['source'] == 'diff_exp':
                source = 'diff_exp'

        # 两个geneset_id到sg_task返回的project_sn等信息是一样的
        task_info = self.ref_rna_v2.get_task_info(geneset_info['task_id'])

        table_name = "KEGG"
        collection_name = "sg_geneset_kegg_class"
        # ref_rna_v2_new.geneset
        to_file = ['ref_rna_v2_new.geneset.export_multi_gene_list(geneset_kegg)',
                   "ref_rna_v2.export_kegg_table2(kegg_table)",
                   "ref_rna_v2.export_kegg_level_table(kegg_table_2)",
                   "ref_rna_v2.export_add_info(add_info)"]

        # geneset_type只有2种，T或者G,s而页面限制只能同时取G或者T，所以这里用一个geneset_id作为示例即可
        option = {
                "task_id": data.task_id,
                "geneset_kegg": json.dumps({'geneset_id': data.geneset_id, 'source': source}),
                "kegg_table": data.geneset_id.split(",")[0],
                "kegg_table_2": data.geneset_id.split(",")[0],
                "add_info": geneset_info['task_id'] + "\t" + data.geneset_type}

        main_table_name = table_name + "_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S")

        mongo_data = [
            ('project_sn', task_info['project_sn']),
            ('task_id', task_info['task_id']),
            ('status', 'start'),
            ('name', main_table_name),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ("version","v3"),
            ("params", json.dumps(params_json, sort_keys=True, separators=(',', ':'))),
        ]
        main_table_id = self.ref_rna_v2.insert_main_table(collection_name, mongo_data)
        new_task_id = self.ref_rna_v2.get_new_id(task_info['task_id'])
        main_table_data = {'run_id': new_task_id}
        update_info = {str(main_table_id): collection_name}

        options = {
            "submit_location": data.submit_location,
            "task_type": data.task_type,
            'main_table_data': main_table_data,
            'update_info': json.dumps(update_info),
            "main_table_id": str(main_table_id),
            "geneset_id": data.geneset_id,
            "geneset_type": data.geneset_type,
            "type": data.type,
            "source": source
        }
        options.update(option)
        if "database_version" in  task_info:
            kegg_version = task_info["database_version"].get("kegg", "")
            options.update({"kegg_version": kegg_version})
        # print(options)
        self.set_sheet_data(name=task_name, options=options, main_table_name=main_table_name, module_type="workflow",
                            to_file=to_file, project_sn=task_info['project_sn'], task_id=task_info['task_id'],
                            new_task_id=new_task_id)

        task_info = super(GenesetKeggAction, self).POST()

        task_info['content'] = {
            'ids': {
                'id': str(main_table_id),
                'name': main_table_name
            }}
        geneset_info = self.ref_rna_v2.insert_geneset_info(data.geneset_id, collection_name, str(main_table_id))
        if 'group_id' in data and str(data.group_id).lower() != 'all':
            _ = self.ref_rna_v2.update_group_is_use(data.task_id, data.group_id)
        if 'control_id' in data:
            _ = self.ref_rna_v2.update_group_compare_is_use(data.task_id, data.control_id)
        return json.dumps(task_info)


if __name__ == '__main__':
    class TestFunction(unittest.TestCase):
        """
        This is test for the tool. Just run this script to do test.
        """

        def test_this(self):
            cmd = 'python /mnt/ilustre/users/sanger-dev/biocluster/bin/webapitest.py '
            cmd += 'post '
            cmd += "-fr no "
            cmd += '-c {} '.format("client03")
            cmd += "s/ref_rna_v2/geneset_kegg "
            cmd += "-b http://192.168.12.101:9090 "
            args = dict(
                submit_location="genesetkegg",
                db_type="ref_rna_v2",
                geneset_id="5cc694d217b2bf5da04d69c3",
                geneset_type="G",
                anno_type="kegg",
                # geneset_kegg="5b0bfa56a4e1af2d43a6c12a",
                task_id="tsg_34011",
                task_type="2",
                type="origin"

            )

            arg_names, arg_values = args.keys(), args.values()
            cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
            print(cmd)
            os.system(cmd)


    unittest.main()

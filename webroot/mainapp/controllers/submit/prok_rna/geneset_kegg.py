# -*- coding: utf-8 -*-
# __author__ = 'qindanhua'
import web
import json
import datetime
import unittest
from mainapp.controllers.project.prok_rna_controller import ProkRNAController
from mbio.api.to_file.prok_rna import *
from mainapp.libs.signature import check_sig


class GenesetKeggAction(ProkRNAController):
    def __init__(self):
        super(GenesetKeggAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        default_argu = ['geneset_id', 'submit_location', 'type', 'task_type','task_id']
        for argu in default_argu:
            if not hasattr(data, argu):
                info = {'success': False, 'info': '%s参数缺少!' % argu}
                return json.dumps(info)
        if hasattr(data, 'submit_type'):
            submit_type = int(data.submit_type)
        else:
            submit_type = 0

        task_name = 'prok_rna.report.geneset_kegg'
        params_json = {
            "submit_location": data.submit_location,
            "task_type": int(data.task_type),
            "geneset_id": data.geneset_id,
            "task_id": data.task_id,
            "type": data.type,
        }
        # 判断传入的基因集id是否存在
        geneset_info = {}
        geneset_documents = []

        for geneset in data.geneset_id.split(","):
            geneset_info = self.prok_rna.get_main_info(geneset, 'sg_geneset', data.task_id)
            if not geneset_info:
                info = {"success": False, "info": "geneset不存在，请确认参数是否正确！!"}
                return json.dumps(info)
            geneset_documents.append(geneset_info)

        source = 'non_diff_exp'
        if len(geneset_documents) == 1:
            my_result = geneset_documents[0]
            if 'source' in my_result and my_result['source'] == 'diff_exp':
                source = 'diff_exp'
        # 两个geneset_id到sg_task返回的project_sn等信息是一样的
        task_info = self.prok_rna.get_task_info(geneset_info['task_id'])

        table_name = "Kegg"
        collection_name = "sg_geneset_kegg_class"
        to_file = ['prok_rna.export_multi_gene_list(geneset_kegg)',
                   "prok_rna.export_kegg_table(kegg_table)",
                   "prok_rna.export_kegg_level_table(kegg_table_2)",
                   "prok_rna.export_add_info(add_info)"]

        option = {"geneset_kegg": data.geneset_id,
                  "kegg_table": data.geneset_id.split(",")[0],
                  "kegg_table_2": data.geneset_id.split(",")[0],
                  "add_info":geneset_info['task_id']}

        main_table_name = 'Geneset' + table_name + "Class_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S")

        mongo_data = [
            ('project_sn', task_info['project_sn']),
            ('task_id', task_info['task_id']),
            ('status', 'start'),
            ('name', main_table_name),
            ('submit_type', submit_type),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ("params", json.dumps(params_json, sort_keys=True, separators=(',', ':'))),
        ]
        main_table_id = self.prok_rna.insert_main_table(collection_name, mongo_data)
        update_info = {str(main_table_id): collection_name}

        options = {
            "task_id": data.task_id,
            "submit_location": data.submit_location,
            "task_type": data.task_type,
            'update_info': json.dumps(update_info),
            "main_table_id": str(main_table_id),
            "geneset_id": data.geneset_id,
            "type":data.type,
            "source": source
        }
        options.update(option)

        if "database_version" in task_info:
            kegg_version = task_info["database_version"].get("kegg", "2017")
            options.update({"kegg_version": kegg_version})
        # print(options)
        self.set_sheet_data(name=task_name, options=options, main_table_name=main_table_name, module_type="workflow",
                            to_file=to_file, project_sn=task_info['project_sn'], task_id=task_info['task_id'])

        task_info = super(GenesetKeggAction, self).POST()

        task_info['content'] = {
            'ids': {
                'id': str(main_table_id),
                'name': main_table_name
            }
        }
        geneset_info = self.prok_rna.insert_geneset_info(data.geneset_id, collection_name, str(main_table_id))
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
        cmd += "s/prok_rna/geneset_kegg "
        cmd += "-b http://192.168.12.102:9090 "
        args = dict(
            task_id="prok_rna_srna",
            task_type="2",
            submit_location="genesetkegg",
            geneset_id="5b70dcb5a4e1af118c63916c,5b70dd58a4e1af13f89b4362",
            type="latest"
        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)


if __name__ == '__main__':
    unittest.main()

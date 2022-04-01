# -*- coding: utf-8 -*-
# __author__ = 'qindanhua'
import web
import json
import datetime
from mainapp.controllers.project.smallrna_controller import SmallRnaController
from mbio.api.to_file.smallrna import *
from mainapp.libs.signature import check_sig
import os
import unittest
from biocluster.config import Config

class GenesetKeggAction(SmallRnaController):
    def __init__(self):
        super(GenesetKeggAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        default_argu = ['geneset_id', 'geneset_type', 'submit_location', 'type', 'task_type', 'task_id']
        for argu in default_argu:
            if not hasattr(data, argu):
                info = {'success': False, 'info': '%s参数缺少!' % argu}
                return json.dumps(info)

        task_name = 'smallrna.report.geneset_kegg'
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
        for geneset in data.geneset_id.split(","):
            geneset_info = self.smallrna.get_main_info(geneset, 'sg_geneset', data.task_id)
            if not geneset_info:
                info = {"success": False, "info": "geneset不存在，请确认参数是否正确！!"}
                return json.dumps(info)
        # 两个geneset_id到sg_task返回的project_sn等信息是一样的
        task_info = self.smallrna.get_task_info(geneset_info['task_id'])

        table_name = "Kegg"
        collection_name = "sg_geneset_kegg_class"
        to_file = ['smallrna.export_multi_gene_list(geneset_kegg)',
                   "smallrna.export_kegg_table(kegg_table)",
                   "smallrna.export_kegg_level_table(kegg_table_2)",
                   "smallrna.export_add_info(add_info)"]

        # geneset_type只有2种，T或者G,s而页面限制只能同时取G或者T，所以这里用一个geneset_id作为示例即可
        option = {"geneset_kegg": data.geneset_id,
                  "kegg_table": data.geneset_id.split(",")[0],
                  "kegg_table_2": data.geneset_id.split(",")[0],
                  "add_info": geneset_info['task_id'] + "\t" + data.geneset_type}

        main_table_name = 'Geneset' + table_name + "Class_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]

        mongo_data = [
            ('project_sn', task_info['project_sn']),
            ('task_id', task_info['task_id']),
            ('status', 'start'),
            ('name', main_table_name),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ("params", json.dumps(params_json, sort_keys=True, separators=(',', ':'))),
        ]
        main_table_id = self.smallrna.insert_main_table(collection_name, mongo_data)
        update_info = {str(main_table_id): collection_name}

        options = {
            "submit_location": data.submit_location,
            "task_type": data.task_type,
            'update_info': json.dumps(update_info),
            "main_table_id": str(main_table_id),
            "geneset_id": data.geneset_id,
            "geneset_type": data.geneset_type,
            "type": data.type,
        }
        options.update(option)
        if "database_version" in task_info:
            kegg_version = task_info["database_version"].get("kegg", "")
            options.update({"kegg_version": kegg_version})
        elif "genome_id" in task_info:
            # 判断基因组注释版本
            database = Config().get_mongo_client(mtype='ref_rna_v2', dydb_forbid=True)[Config().get_mongo_dbname('ref_rna_v2', dydb_forbid=True)]
            collection = database['sg_genome_db']
            genome_doc = collection.find_one({'genome_id': task_info["genome_id"]})
            if "database_version" in genome_doc and "kegg" in genome_doc["database_version"]:
                options.update({
                    "kegg_version": genome_doc["database_version"]["kegg"]
                })
        # print(options)
        self.set_sheet_data(name=task_name, options=options, main_table_name=main_table_name, module_type="workflow",
                            to_file=to_file, project_sn=task_info['project_sn'], task_id=task_info['task_id'])

        task_info = super(GenesetKeggAction, self).POST()

        task_info['content'] = {
            'ids': {
                'id': str(main_table_id),
                'name': main_table_name
            }}
        geneset_info = self.smallrna.insert_geneset_info(data.geneset_id, collection_name, str(main_table_id))
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
        cmd += "s/smallrna/geneset_kegg "
        cmd += "-b http://192.168.12.102:9090 "
        args = dict(
            submit_location="genesetkegg",
            db_type="smallrna",
            geneset_id="5b0686eca4e1af2636e97f7e,5b06899da4e1af3a72687201",
            geneset_type="T",
            anno_type="kegg",
            geneset_kegg="5b0bfa56a4e1af2d43a6c12a",
            task_id="RefrnaV2_7320",
            task_type="2",
            type="origin"

        )

        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)


if __name__ == '__main__':
    unittest.main()

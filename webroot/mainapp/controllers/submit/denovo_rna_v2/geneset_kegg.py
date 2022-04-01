# -*- coding: utf-8 -*-
# __author__ = 'qindanhua'
import web
import json
import datetime
from mainapp.controllers.project.denovo_rna_v2_controller import DenovoRnaV2Controller
from mbio.api.to_file.denovo_rna_v2 import *
from mainapp.libs.signature import check_sig
import unittest


class GenesetKeggAction(DenovoRnaV2Controller):
    def __init__(self):
        super(GenesetKeggAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        default_argu = ['geneset_id', 'geneset_type', 'submit_location', 'type', 'task_type','task_id']
        for argu in default_argu:
            if not hasattr(data, argu):
                var = []
                var.append(argu)
                info = {'success': False, 'info': "Lack argument: %s" % (argu), "code": 'C1601401', "variables": var}
                return json.dumps(info)


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
            geneset_info = self.denovo_rna_v2.get_main_info(geneset, 'sg_geneset', data.task_id)
            if not geneset_info:
                info = {"success": False, "info": "geneset不存在，请确认参数是否正确！!", "code": 'C1601401', "variables": ''}
                return json.dumps(info)
        # 两个geneset_id到sg_task返回的project_sn等信息是一样的
        task_info = self.denovo_rna_v2.get_task_info(geneset_info['task_id'])


        if "version" in task_info:
             version=task_info["version"]
        else:
            version="v1"
        if version == "v1":
            task_name = 'denovo_rna_v2.report.geneset_kegg'
            table_name = "Kegg"
            collection_name = "sg_geneset_kegg_class"
            to_file = ['denovo_rna_v2.export_multi_gene_list(geneset_kegg)',
                       "denovo_rna_v2.export_kegg_table(kegg_table)",
                       "denovo_rna_v2.export_kegg_level_table(kegg_table_2)",
                       "denovo_rna_v2.export_add_info(add_info)"]

            # geneset_type只有2种，T或者G,s而页面限制只能同时取G或者T，所以这里用一个geneset_id作为示例即可
            option = {"geneset_kegg": data.geneset_id,
                      "kegg_table": data.geneset_id.split(",")[0],
                      "kegg_table_2": data.geneset_id.split(",")[0],
                      "add_info":geneset_info['task_id'] + "\t" + data.geneset_type}

            main_table_name = 'Geneset' + table_name + "Class_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S")

            mongo_data = [
                ('project_sn', task_info['project_sn']),
                ('anno_detail', "v2"),
                ('task_id', task_info['task_id']),
                ('status', 'start'),
                ('name', main_table_name),
                ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
                ("params", json.dumps(params_json, sort_keys=True, separators=(',', ':'))),
            ]
            main_table_id = self.denovo_rna_v2.insert_main_table(collection_name, mongo_data)
            update_info = {str(main_table_id): collection_name}

            options = {
                "task_id": data.task_id,
                "submit_location": data.submit_location,
                "task_type": data.task_type,
                'update_info': json.dumps(update_info),
                "main_table_id": str(main_table_id),
                "geneset_id": data.geneset_id,
                "geneset_type": data.geneset_type,
                "type":data.type,
                }
            options.update(option)

        else:
            task_name = 'denovo_rna_v3.report.geneset_kegg'
            table_name = "Kegg"
            collection_name = "sg_geneset_kegg_class"
            if 'source' in geneset_info and geneset_info['source'] == 'diff_exp':
                source = 'diff_exp'
            else:
                source = 'non_diff_exp'
            to_file = ['denovo_rna_v3.export_multi_gene_list(geneset_kegg)',
                       "denovo_rna_v3.export_kegg_table(kegg_table)",
                       "denovo_rna_v3.export_kegg_level_table(kegg_table_2)",
                       "denovo_rna_v3.export_add_info(add_info)"]

            # geneset_type只有2种，T或者G,s而页面限制只能同时取G或者T，所以这里用一个geneset_id作为示例即可
            option = {"geneset_kegg": json.dumps({'geneset_id': data.geneset_id, 'source': source}),
                      "kegg_table": data.geneset_id.split(",")[0],
                      "kegg_table_2": data.geneset_id.split(",")[0],
                      "add_info": geneset_info['task_id'] + "\t" + data.geneset_type,
                      'source':source
                         }

            main_table_name = 'Geneset' + table_name + "Class_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S")

            mongo_data = [
                ('project_sn', task_info['project_sn']),
                ('version', 'v2'),
                ('anno_detail', "v2"),
                ('task_id', task_info['task_id']),
                ('status', 'start'),
                ('name', main_table_name),
                ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
                ("params", json.dumps(params_json, sort_keys=True, separators=(',', ':'))),
            ]
            main_table_id = self.denovo_rna_v2.insert_main_table(collection_name, mongo_data)
            update_info = {str(main_table_id): collection_name}

            options = {
                "task_id": data.task_id,
                "submit_location": data.submit_location,
                "task_type": data.task_type,
                'update_info': json.dumps(update_info),
                "main_table_id": str(main_table_id),
                "geneset_id": data.geneset_id,
                "geneset_type": data.geneset_type,
                "type": data.type,
            }
            options.update(option)

            if "database_version" in  task_info:
                kegg_version = task_info["database_version"].get("kegg", "")
                options.update({"kegg_version": kegg_version})
        # print(options)
        self.set_sheet_data(name=task_name, options=options, main_table_name=main_table_name, module_type="workflow",
                            to_file=to_file, project_sn=task_info['project_sn'], task_id=task_info['task_id'])

        task_info = super(GenesetKeggAction, self).POST()

        task_info['content'] = {
            'ids': {
                'id': str(main_table_id),
                'name': main_table_name
                }}
        geneset_info = self.denovo_rna_v2.insert_geneset_info(data.geneset_id, collection_name, str(main_table_id))
        if 'group_id' in data and str(data.group_id).lower() != 'all':
            _ = self.denovo_rna_v2.update_group_is_use(data.task_id, data.group_id)
        if 'control_id' in data:
            _ = self.denovo_rna_v2.update_group_compare_is_use(data.task_id, data.control_id)
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
            cmd += "s/denovo_rna_v2/geneset_kegg "
            cmd += "-b http://bcl.tsg.com "
            args = dict(
                submit_location="genesetkegg",
                db_type="denovo_rna_v2",
                # geneset_id="5d7992c517b2bf1d40fdb8c5",#5c7e9f0417b2bf7eb4e918dd",  # ObjectId("5cb952d28876b737418b4568")
                geneset_id="5d7992b717b2bf1d40fd0be3",
                geneset_type="T",
                anno_type="kegg",
                # geneset_kegg="5b0bfa56a4e1af2d43a6c12a",
                task_id="tsg_35544",
                task_type="2",
                type="origin"
            )
            arg_names, arg_values = args.keys(), args.values()
            cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
            print(cmd)
            os.system(cmd)


    unittest.main()

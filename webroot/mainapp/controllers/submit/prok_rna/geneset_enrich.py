# -*- coding: utf-8 -*-
# __author__ = 'qindanhua'
# -*- coding: utf-8 -*-
import web
import json
import datetime
import unittest
from collections import OrderedDict
from mainapp.controllers.project.prok_rna_controller import ProkRNAController
from mbio.api.to_file.prok_rna import *
from mainapp.libs.signature import check_sig


class GenesetEnrichAction(ProkRNAController):
    def __init__(self):
        super(GenesetEnrichAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        default_argu = ['geneset_id', 'submit_location', 'task_type', 'anno_type', 'method', "type"]
        for argu in default_argu:
            if not hasattr(data, argu):
                info = {'success': False, 'info': '%s参数缺少!' % argu}
                return json.dumps(info)
        if hasattr(data, 'submit_type'):
            submit_type = int(data.submit_type)
        else:
            submit_type = 0

        task_name = 'prok_rna.report.geneset_enrich'
        task_type = 'workflow'

        params_json = {
            "submit_location": data.submit_location,
            "task_type": int(data.task_type),
            "anno_type": data.anno_type,
            "method": data.method,
            "type": data.type,
            "geneset_id": data.geneset_id,
            "task_id": data.task_id,
        }
        # 判断传入的基因集id是否存在
        geneset_info = self.prok_rna.get_main_info(data.geneset_id, 'sg_geneset', data.task_id)
        if not geneset_info:
            info = {"success": False, "info": "geneset不存在，请确认参数是否正确！!"}
            return json.dumps(info)
        task_info = self.prok_rna.get_task_info(geneset_info['task_id'])
        geneset_name = geneset_info['name']

        if 'source' in geneset_info and geneset_info['source'] == 'diff_exp':
            source = 'diff_exp'
        else:
            source = 'non_diff_exp'
        if data.anno_type == "go":
            table_name = "Go"
            collection_name = "sg_geneset_go_enrich"
            to_file = ['prok_rna.export_gene_list(geneset_list)',
                       # 'prok_rna.export_all_list(all_list)',
                       'prok_rna.export_go_list(go_list)']
            infile = {"go_list": data.geneset_id}
        elif data.anno_type == "kegg":
            table_name = "Kegg"
            collection_name = "sg_geneset_kegg_enrich"
            to_file = ['prok_rna.export_gene_list(geneset_list)',
                       'prok_rna.export_multi_gene_list(geneset_kegg)',
                       'prok_rna.export_all_list(all_list)',
                       'prok_rna.export_kegg_table(kegg_table)',
                       "prok_rna.export_add_info(add_info)"]
            infile = {"kegg_table": data.geneset_id, "add_info": geneset_info['task_id']}
        else:
            info = {'success': False, 'info': '不支持该富集分子!'}
            return json.dumps(info)

        main_table_name = 'Geneset' + table_name + "Enrich_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S")

        mongo_data = [
            ('project_sn', task_info['project_sn']),
            ('task_id', task_info['task_id']),
            ('status', 'start'),
            ('name', main_table_name),
            ('submit_type', submit_type),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ("params", json.dumps(params_json, sort_keys=True, separators=(',', ':')))
        ]
        main_table_id = self.prok_rna.insert_main_table(collection_name, mongo_data)
        update_info = {str(main_table_id): collection_name}

        options = {
            "task_id": data.task_id,
            "geneset_kegg": data.geneset_id,
            "submit_location": data.submit_location,
            "task_type": data.task_type,
            'update_info': json.dumps(update_info),
            "main_table_id": str(main_table_id),
            "anno_type": data.anno_type,
            "method": data.method,
            "geneset_list": data.geneset_id,
            "type": data.type,
            "geneset_id": data.geneset_id,
            "all_list": data.geneset_id,
            "source": source,
            "geneset_name": geneset_name,
        }
        options.update(infile)
        if "database_version" in  task_info:
            kegg_version = task_info["database_version"].get("kegg", "")
            options.update({"kegg_version": kegg_version})

        if "database_version" in task_info and data.anno_type == "kegg":
            kegg_version = task_info["database_version"].get("kegg", "2017")
            options.update({"kegg_version": kegg_version})
        if "database_version" in task_info and data.anno_type == "go":
            go_version = task_info["database_version"].get("go", "20200628").split("_")[0]
            options.update({"go_version": go_version})

        self.set_sheet_data(name=task_name, options=options, main_table_name=main_table_name, module_type=task_type,
                            to_file=to_file, task_id=task_info["task_id"], project_sn=task_info["project_sn"])
        task_info = super(GenesetEnrichAction, self).POST()
        task_info['content'] = {
            'ids': {
                'id': str(main_table_id),
                'name': main_table_name
            }
        }
        self.prok_rna.insert_geneset_info(data.geneset_id, collection_name, str(main_table_id))
        if 'group_id' in data and str(data.group_id).lower() != 'all':
            _ = self.prok_rna.update_group_is_use(data.task_id, data.group_id)
        if 'control_id' in data:
            _ = self.prok_rna.update_group_compare_is_use(data.task_id, data.control_id)

        return json.dumps(task_info)

class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    import os
    def test_this(self):
        cmd = 'python /mnt/ilustre/users/sanger-dev/biocluster/bin/webapitest.py '
        cmd += 'post '
        cmd += "-fr no "
        cmd += '-c {} '.format("client03")
        cmd += "s/prok_rna/geneset_enrich "
        cmd += "-b http://192.168.12.102:9090 "
        args = dict(
            task_id="prok_rna",
            task_type="2",
            submit_location="geneset_keggenrich",
            geneset_id="5b70dd58a4e1af13f89b4362",
            type="origin",
            anno_type="kegg",
            method='BH'
        )
        '''
        args = dict(
            task_id="prok_rna",
            task_type="2",
            submit_location="genesetgo",
            geneset_id="5a9de1715b8c915efb2a3293",
            type="origin",
            anno_type="go"
        )
        '''
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)


if __name__ == '__main__':
    unittest.main()


# -*- coding: utf-8 -*-
# __author__ = 'qindanhua'

import datetime
import json

import web

from mainapp.controllers.project.ref_rna_controller import RefRnaController
from mainapp.libs.signature import check_sig


class GenesetEnrichAction(RefRnaController):
    def __init__(self):
        super(GenesetEnrichAction, self).__init__()

    @check_sig
    def POST(self):
        data = web.input()
        default_argu = ['geneset_id', 'geneset_type', 'submit_location', 'anno_type', 'method']
        for argu in default_argu:
            if not hasattr(data, argu):
                info = {'success': False, 'info': '%s参数缺少!' % argu}
                return json.dumps(info)

        task_name = 'ref_rna.report.geneset_enrich'
        task_type = 'workflow'

        params_json = {
            "submit_location": data.submit_location,
            "task_type": data.task_type,
            "geneset_id": data.geneset_id,
            "anno_type": data.anno_type,
            "method": data.method,
            "geneset_type": data.geneset_type
        }
        geneset_info = self.ref_rna.get_main_info(data.geneset_id, 'sg_geneset')
        if not geneset_info:
            info = {"success": False, "info": "geneset不存在，请确认参数是否正确！!"}
            return json.dumps(info)
        task_info = self.ref_rna.get_task_info(geneset_info['task_id'])

        if data.anno_type == "go":
            table_name = "Go"
            collection_name = "sg_geneset_go_enrich"
            to_file = ['ref_rna.export_gene_list(genset_list)', 'ref_rna.export_all_list(all_list)',
                       'ref_rna.export_go_list(go_list)']
            infile = {"go_list": data.geneset_id}
        elif data.anno_type == "kegg":
            table_name = "Kegg"
            collection_name = "sg_geneset_kegg_enrich"
            to_file = ['ref_rna.export_gene_list(genset_list)', 'ref_rna.export_all_list(all_list)',
                       'ref_rna.export_kegg_table(kegg_table)', "ref_rna.export_add_info(add_info)"]
            infile = {"kegg_table": data.geneset_id, "add_info": geneset_info['task_id'] + "\t" + data.geneset_type}
        else:
            info = {'success': False, 'info': '不支持该富集分子!'}
            return json.dumps(info)

        main_table_name = 'Geneset' + table_name + "Enrich" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]

        mongo_data = [
            ('project_sn', task_info['project_sn']),
            ('task_id', task_info['task_id']),
            ('status', 'start'),
            ('name', main_table_name),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ("params", json.dumps(params_json, sort_keys=True, separators=(',', ':')))
        ]
        main_table_id = self.ref_rna.insert_main_table(collection_name, mongo_data)
        update_info = {str(main_table_id): collection_name}

        options = {
            "submit_location": data.submit_location,
            "task_type": data.task_type,
            'update_info': json.dumps(update_info),
            "main_table_id": str(main_table_id),
            "geneset_type": data.geneset_type,
            "anno_type": data.anno_type,
            "method": data.method,
            "genset_list": data.geneset_id,
            "all_list": data.geneset_id
        }
        options.update(infile)

        self.set_sheet_data(name=task_name, options=options, main_table_name=main_table_name, module_type=task_type,
                            to_file=to_file, task_id=task_info["task_id"], project_sn=task_info["project_sn"])
        task_info = super(GenesetEnrichAction, self).POST()
        task_info['content'] = {'ids': {'id': str(main_table_id), 'name': main_table_name}}

        return json.dumps(task_info)

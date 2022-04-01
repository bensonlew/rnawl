# -*- coding: utf-8 -*-
# __author__ = 'qindanhua'
import web
import json
import datetime
from mainapp.controllers.project.ref_rna_controller import RefRnaController
from mbio.api.to_file.ref_rna import *
from mainapp.libs.signature import check_sig


class GenesetKeggAction(RefRnaController):
    def __init__(self):
        super(GenesetKeggAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        default_argu = ['geneset_id', 'geneset_type', 'submit_location']
        for argu in default_argu:
            if not hasattr(data, argu):
                info = {'success': False, 'info': '%s参数缺少!' % argu}
                return json.dumps(info)

        task_name = 'ref_rna.report.geneset_kegg'
        task_type = ''
        params_json = {
            "submit_location": data.submit_location,
            "task_type": data.task_type,
            "geneset_id": data.geneset_id,
            "geneset_type": data.geneset_type
        }
        # 判断传入的基因集id是否存在
        geneset_info = {}
        for geneset in data.geneset_id.split(","):
            geneset_info = self.ref_rna.get_main_info(geneset, 'sg_geneset')
            if not geneset_info:
                info = {"success": False, "info": "geneset不存在，请确认参数是否正确！!"}
                return json.dumps(info)
        task_info = self.ref_rna.get_task_info(geneset_info['task_id'])

        table_name = "Kegg"
        collection_name = "sg_geneset_kegg_class"
        to_file = ['ref_rna.export_multi_gene_list(geneset_kegg)', "ref_rna.export_kegg_table(kegg_table)",
                   # "ref_rna.export_kegg_pdf(kegg_pics)",
                   "ref_rna.export_add_info(add_info)"]
        option = {"geneset_kegg": data.geneset_id, "kegg_table": data.geneset_id.split(",")[0],
                  # "kegg_pics": geneset_info['task_id'] + "\t" + data.geneset_type,
                  "add_info":geneset_info['task_id'] + "\t" + data.geneset_type}

        main_table_name = 'Geneset' + table_name + "Class_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]

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
            "task_type": "workflow",
            'update_info': json.dumps(update_info),
            "main_table_id": str(main_table_id),
            "geneset_id": data.geneset_id,
            "geneset_type": data.geneset_type
            }
        options.update(option)
        # print(options)
        self.set_sheet_data(name=task_name, options=options, main_table_name=main_table_name, module_type="workflow",
                            to_file=to_file, project_sn=task_info['project_sn'], task_id=task_info['task_id'])

        task_info = super(GenesetKeggAction, self).POST()

        task_info['content'] = {
            'ids': {
                'id': str(main_table_id),
                'name': main_table_name
                }}
        geneset_info = self.ref_rna.insert_geneset_info(data.geneset_id, collection_name, str(main_table_id))
        #if geneset_info:
            #print "geneset_info插入成功"
        return json.dumps(task_info)

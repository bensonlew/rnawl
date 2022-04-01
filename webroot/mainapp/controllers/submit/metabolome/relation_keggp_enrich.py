# -*- coding: utf-8 -*-
# __author__ = 'zhaoyuzhou'
# last_modifiy = modified 2021.11.30

import web
import json
import datetime
from mainapp.controllers.project.metabolome_controller import MetabolomeController
from mainapp.models.mongo.metabolome import Metabolome
from mainapp.libs.signature import check_sig


class RelationKeggpEnrichAction(MetabolomeController):
    def __init__(self):
        super(RelationKeggpEnrichAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        print "data:", data
        default_argu = ["submit_location", 'task_type', 'task_id', 'metab_keggp_table', 'metab_set_table', 'trans_keggp_main_id',
                        'trans_geneset_main_id', 'padjust_method']
        for arg in default_argu:
            if not hasattr(data, arg):
                info = {'success': False, 'info': 'Argument:{} required'.format(arg)}
                return json.dumps(info)
        metabolome = Metabolome()
        task_name = 'metabolome.report.relation_keggp_enrich'
        module_type = 'workflow'
        task_id = data.task_id
        project_sn, project_type, member_id = metabolome.get_project_info(task_id)
        name = "RelationKeggpEnrich_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        if project_type == "LC":
            metab_desc_path = metabolome.get_metab_desc1(task_id, 'mix')
        elif project_type == "GC":
            metab_desc_path = metabolome.get_metab_desc1(task_id)
        params_json = {
            'metab_keggp_table': data.metab_keggp_table,
            'metab_set_table': data.metab_set_table,
            'trans_keggp_main_id': data.trans_keggp_main_id,
            'trans_geneset_main_id': data.trans_geneset_main_id,
            'padjust_method':data.padjust_method,
            'task_type': int(data.task_type),
            'submit_location': data.submit_location,
            'task_id': task_id
        }
        mongo_data = [
            ('project_sn', project_sn),
            ('task_id', task_id),
            ('status', 'start'),
            ("name", name),
            ("desc", "Running"),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        ]
        mongo_data.append(('params', json.dumps(params_json, sort_keys=True, separators=(',', ':'))))
        main_table_id = metabolome.insert_main_table('relation_keggp_enrich', mongo_data)
        metabolome.insert_main_id('relation_keggp_enrich', main_table_id)
        options = {
            'metab_set_table': data.metab_set_table,
            'anno_overview': data.task_id,
            'ko_overview': data.task_id,
            "metab_desc": metab_desc_path,
            'trans_keggp_main_id': data.trans_keggp_main_id,
            'trans_geneset_main_id': data.trans_geneset_main_id,
            'padjust_method': data.padjust_method,
            'main_table_id': str(main_table_id),
            "name": name
        }
        species = metabolome.get_kegg_species(data.task_id)
        if species:
            options["species"] = species
        else:
            info = {'success': False, 'info': 'not result'}
            return json.dumps(info)
        update_info = {str(main_table_id): 'relation_keggp_enrich'}
        options["update_info"] = json.dumps(update_info)
        to_file = []
        to_file.append('metabolome.export_metab_set1(metab_set_table)')
        to_file.append('metabolome.export_overview(anno_overview)')
        to_file.append('metabolome.export_overview_ko(ko_overview)')
        m_table_name = "Relation/"+name.strip().split("_")[0] + '/' + name
        self.set_sheet_data(name=task_name, options=options, main_table_name=m_table_name,
                            module_type=module_type, project_sn=project_sn, to_file=to_file,
                            task_id=task_id, params=params_json)
        task_info = super(RelationKeggpEnrichAction, self).POST()
        if task_info['success']:
            task_info['content'] = {'ids': {'id': str(main_table_id), 'name': name}}
        return json.dumps(task_info)

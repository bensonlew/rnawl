# -*- coding: utf-8 -*-
# __author__ = 'zhaoyuzhou'
# last_modifiy = modified 2021.11.30

import web
import json
import datetime
from mainapp.controllers.project.metabolome_controller import MetabolomeController
from mainapp.models.mongo.metabolome import Metabolome
from mainapp.libs.signature import check_sig


class RelationO2plsAction(MetabolomeController):
    """ O2pls.controller"""
    def __init__(self):
        super(RelationO2plsAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        global metab_table_path
        data = web.input()
        print "data:", data
        default_argu = ["submit_location", 'task_type', 'task_id', 'metab_table',
                        'trans_exp_main_id', 'group_detail', "scale", "sort_name"]
        for arg in default_argu:
            if not hasattr(data, arg):
                info = {"success": False, "info": "PARAMETERS MISSING: %s" % arg}
                return json.dumps(info)
        metabolome = Metabolome()
        task_name = 'metabolome.report.relation_o2pls'
        module_type = 'workflow'
        task_id = data.task_id
        project_sn, project_type, table_name = metabolome.get_metab_info(data.metab_table, task_id)
        name = "RelationO2pls_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        if not data.scale in ["UV", "Ctr", "Par", "none"]:
            info = {'success': False, 'info': 'PARAMETERS ERROR: wrong value of scale'}
            return json.dumps(info)
        if project_type == "LC":
            metab_table_path = metabolome.get_metab_table(data.metab_table, task_id, 'mix')
            metab_desc_path = metabolome.get_metab_desc(data.metab_table, task_id, 'mix')
        elif project_type == "GC":
            metab_table_path = metabolome.get_metab_table(data.metab_table, task_id)
            metab_desc_path = metabolome.get_metab_desc(data.metab_table, task_id)

        params_json = {
            'metab_table': data.metab_table,
            'trans_exp_main_id': data.trans_exp_main_id,
            "group_detail": json.loads(data.group_detail),
            'group_id': data.group_id,
            'scale': data.scale,
            'task_type': int(data.task_type),
            'submit_location': data.submit_location,
            'task_id': task_id,
            "sort_name": data.sort_name
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
        main_table_id = metabolome.insert_main_table('relation_o2pls', mongo_data)
        metabolome.insert_main_id('relation_o2pls', main_table_id)
        options = {
            'metab_table': metab_table_path,
            "metab_desc": metab_desc_path,
            'trans_exp_main_id': data.trans_exp_main_id,
            "group_detail": data.group_detail,
            'scale': data.scale,
            'main_table_id': str(main_table_id),
            "name": name,
            "sort_name": data.sort_name
        }
        if table_name == "MetabTable_Origin":
            scale = metabolome.get_scale_type(data.metab_table, task_id)
            if scale and scale == "none":
                options["log10"] = True
            else:
                options["log10"] = False
        elif table_name == "raw":
            options["log10"] = True
        update_info = {str(main_table_id): 'relation_o2pls'}
        options["update_info"] = json.dumps(update_info)
        to_file = []
        to_file.append('metabolome.export_group_by_detail(group_detail)')
        m_table_name = "Relation/"+name.strip().split("_")[0] + '/' + name
        self.set_sheet_data(name=task_name, options=options, main_table_name=m_table_name, module_type=module_type,
                            project_sn=project_sn, to_file=to_file, task_id=task_id, params=params_json)
        task_info = super(RelationO2plsAction, self).POST()
        if task_info['success']:
            task_info['content'] = {'ids': {'id': str(main_table_id), 'name': name}}
        return json.dumps(task_info)
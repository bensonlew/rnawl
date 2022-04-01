# -*- coding: utf-8 -*-
# __author__ = 'zhaoyuzhou'
# last_modifiy = modified 2021.11.30

import web
import json
import datetime
from mainapp.controllers.project.metabolome_controller import MetabolomeController
from mainapp.models.mongo.metabolome import Metabolome
from mainapp.libs.signature import check_sig


class RelationProcrustesAction(MetabolomeController):
    """ Procrustes controller"""
    def __init__(self):
        super(RelationProcrustesAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        print "data:", data
        default_argu = ["submit_location", 'task_type', 'task_id', 'metab_table', 'metab_set_table', 'trans_exp_main_id', 'group_id',
                        'trans_geneset_main_id', 'group_detail', "sort_method"]
        for arg in default_argu:
            if not hasattr(data, arg):
                info = {"success": False, "info": "PARAMETERS MISSING: %s" % arg}
                return json.dumps(info)
        metabolome = Metabolome()
        task_name = 'metabolome.report.relation_procrustes'
        module_type = 'workflow'
        task_id = data.task_id
        project_sn, project_type, table_name = metabolome.get_metab_info(data.metab_table, task_id)
        name = "RelationProcrustes_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        if not data.sort_method in ["pca", "pcoa"]:
            info = {'success': False, 'info': 'PARAMETERS ERROR: wrong value of method, pca or pcoa expected!'}
            return json.dumps(info)

        if project_type == "LC":
            metab_table_path = metabolome.get_metab_table(data.metab_table, task_id, 'mix')
        elif project_type == "GC":
            metab_table_path = metabolome.get_metab_table(data.metab_table, task_id)
        params_json = {
            'metab_table': data.metab_table,
            'metab_set_table': data.metab_set_table,
            'trans_exp_main_id':data.trans_exp_main_id,
            'trans_geneset_main_id':data.trans_geneset_main_id,
            "group_detail": json.loads(data.group_detail),
            'group_id': data.group_id,
            'sort_method':data.sort_method,
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
        main_table_id = metabolome.insert_main_table('relation_procrustes', mongo_data)
        metabolome.insert_main_id('relation_procrustes', main_table_id)
        options = {
            'metab_table': metab_table_path,
            'metab_set_table': data.metab_set_table,
            'trans_exp_main_id': data.trans_exp_main_id,
            'trans_geneset_main_id': data.trans_geneset_main_id,
            "group_detail": data.group_detail,
            'method': data.sort_method,
            'main_table_id': str(main_table_id),
            "name": name
        }
        if table_name == "MetabTable_Origin":
            scale = metabolome.get_scale_type(data.metab_table, task_id)
            if scale and scale == "none":
                options["log10"] = True
            else:
                options["log10"] = False
        elif table_name == "raw":
            options["log10"] = True
        update_info = {str(main_table_id): 'relation_procrustes'}
        options["update_info"] = json.dumps(update_info)
        to_file = []
        if data.metab_set_table != "all":
            to_file.append('metabolome.export_metab_set1(metab_set_table)')
        else:
            pass
        to_file.append('metabolome.export_group_by_detail(group_detail)')
        task_info = metabolome.get_task_info(data.task_id)
        m_table_name = "Relation/"+name.strip().split("_")[0] + '/' + name
        self.set_sheet_data(name=task_name, options=options, main_table_name=m_table_name,
                            module_type=module_type, project_sn=project_sn, to_file=to_file,
                            task_id=task_id, params=params_json)
        task_info = super(RelationProcrustesAction, self).POST()
        if task_info['success']:
            task_info['content'] = {'ids': {'id': str(main_table_id), 'name': name}}
        return json.dumps(task_info)
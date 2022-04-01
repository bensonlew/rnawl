# -*- coding: utf-8 -*-

import os
import web
import json
import datetime
from mainapp.controllers.project.metabolome_controller import MetabolomeController
from mainapp.models.mongo.metabolome import Metabolome
from mainapp.libs.signature import check_sig
from bson import ObjectId
import sqlite3
import unittest



class SamplePlsdaAction(MetabolomeController):
    def __init__(self):
        super(SamplePlsdaAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        print data
        default_argu = ['submit_location', 'task_type', 'task_id', 'group_detail','group_id','metab_table','trans','confidence','replace']
        # check arg
        for arg in default_argu:
            if not hasattr(data, arg):
                info = {"success": False, "info": "PARAMETERS MISSING: %s" % arg}
                return json.dumps(info)
        metabolome = Metabolome()
        task_name = 'metabolome.report.sample_plsda'
        module_type = 'workflow'
        task_id = data.task_id
        #project_sn, project_type, member_id = metabolome.get_project_info(task_id)
        project_sn, project_type, table_name = metabolome.get_metab_info(data.metab_table, task_id)
        task_info = self.Metabolome.get_task_info(data.task_id)
        if "save_pdf" in task_info and int(task_info["save_pdf"]) == 1:
            name = "ExpPLSDA_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
            main_table_name = '2.SampleComp/04.ExpPLSDA/' + name
        else:
            name = "Plsda_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
            main_table_name = name.strip().split("_")[0] + '/' + name

        params_json = {
            'metab_table' : data.metab_table,
            'group_id' : data.group_id,
            'group_detail' : json.loads(data.group_detail),
            'task_type': int(data.task_type),
            'submit_location': data.submit_location,
            'task_id': task_id,
            'trans': data.trans,
            'confidence': float(data.confidence),
            'replace': int(data.replace)
        }

        mongo_data = [
            ('project_sn', project_sn),
            ('task_id', task_id),
            ('status', 'start'),
            ("name", name),
            ("desc", "Running"),
            ("version","3.0"),
            ("metab_table", ObjectId(data.metab_table)),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        ]
        if project_type == "LC":
            if hasattr(data, "table_type"):
                table_type= data.table_type
                params_json['table_type'] = table_type
                # info = {"success": False, "info": "LC项目必须输入metab_tabel阴阳离子类型!", 'code':'C2300602'}
                # return json.dumps(info)
            else:
                table_type = 'mix'
            metab_table_path = metabolome.get_metab_table(data.metab_table, task_id, table_type)
            metab_desc_path = metabolome.get_metab_desc(data.metab_table, task_id, table_type)
        elif project_type == "GC":
            table_type = 'pos'
            metab_table_path = metabolome.get_metab_table(data.metab_table, task_id)
            metab_desc_path = metabolome.get_metab_desc(data.metab_table, task_id)

        mongo_data.append(('params', json.dumps(params_json, sort_keys=True, separators=(',', ':'))))
        metabolome = Metabolome()
        main_table_id = metabolome.insert_main_table('sample_plsda', mongo_data)
        metabolome.insert_main_id('sample_plsda', main_table_id)
        options = {
            "metab_table": metab_table_path,
            "group_detail": data.group_detail,
            'main_table_id': str(main_table_id),
            'group': data.group_id,
            "metab_desc": metab_desc_path,
            "confidence": data.confidence,
            "perm": data.replace,
            "data_trans": data.trans,
            "table_type" : table_type,
            "name": name
        }
        if "save_pdf" in task_info and int(task_info["save_pdf"]) == 1:
            options["save_pdf"] = 1
        update_info = {str(main_table_id): 'sample_plsda'}
        options["update_info"] = json.dumps(update_info)
        to_file = []
        to_file.append('metabolome.export_group_by_detail(group)')
        print '{}'.format(main_table_name)
        self.set_sheet_data(name=task_name, options=options, main_table_name=main_table_name,
                            module_type=module_type, project_sn=project_sn,to_file=to_file,
                            task_id=task_id, params=params_json)
        task_info = super(SamplePlsdaAction, self).POST()
        if task_info['success']:
            task_info['content'] = {'ids': {'id': str(main_table_id), 'name': name}}
        return json.dumps(task_info)

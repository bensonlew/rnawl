# -*- coding: utf-8 -*-
import os
import web
import json
import datetime
from mainapp.controllers.project.metabolome_controller import MetabolomeController
from mainapp.models.mongo.metabolome import Metabolome
from mainapp.libs.signature import check_sig
import sqlite3
import unittest
# __author__ = 'shaohua.yuan'


class AnnoHmdbAction(MetabolomeController):
    def __init__(self):
        super(AnnoHmdbAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        for arg in ['metab_table', 'task_id', 'submit_location', 'task_type']:  # table_type
            if not hasattr(data, arg):
                info = {'success': False, 'info': 'Argument:{} required'.format(arg)}
                return json.dumps(info)
        task_name = 'metabolome.report.anno_hmdb'
        module_type = 'workflow'
        params_json = {
            'metab_table': data.metab_table,
            # 'table_type': data.table_type,
            'task_id': data.task_id,
            'submit_location': data.submit_location,
            'task_type': int(data.task_type)
        }
        if hasattr(data, 'table_type'):
            params_json['table_type'] = data.table_type
        task_info = self.Metabolome.get_task_info(data.task_id)
        mongo_data = [
            # ('project_sn', data.project_sn),
            ('project_sn', task_info['project_sn']),
            ('task_id', data.task_id),
            ('status', 'start'),
            ('desc', '正在计算'),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        ]
        main_table_name = "Anno_Hmdb_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        mongo_data.append(('name', main_table_name))
        mongo_data.append(('params', json.dumps(params_json, sort_keys=True, separators=(',', ':'))))
        metabolome = Metabolome()
        v2 = metabolome.conmon_find_one('sg_task', {"task_id": data.task_id})

        if 'version' in v2  and float(v2['version'])>1:
            if v2['type'] == 'LC':
                table_type = 'mix'
            else:
                table_type = 'pos'
        else:
            if hasattr(data, 'table_type'):
                table_type = data.table_type
            else:
                table_type = 'pos'
        main_table_id = metabolome.insert_main_table('anno_hmdb', mongo_data)
        metabolome.insert_main_id('anno_hmdb', main_table_id)
        update_info = {str(main_table_id): 'anno_hmdb'}
        metabtable =  metabolome.get_metab_table(data.metab_table, data.task_id, table_type)
        hmdb_overview = metabolome.get_hmdb_anno(data.task_id)
        if hmdb_overview == "":
            info = {"success": False, "info": "this project has no origin hmdb annotation result" }
            return json.dumps(info)
        options = {
            'update_info': json.dumps(update_info),
            'anno_overview': hmdb_overview,
            "metabtable": metabtable,
            'main_table_id': str(main_table_id),
            'name': main_table_name,
            "type": "annohmdb"
        }
        if "save_pdf" in task_info and int(task_info["save_pdf"]) == 1:
            options["save_pdf"] = 1
            m_table_name = '3.Anno/03.AnnoHmdb/' + main_table_name
        else:
            m_table_name = main_table_name.strip().split("_")[0] + '/' + main_table_name
        self.set_sheet_data(name=task_name, options=options,
                            main_table_name=m_table_name,
                            task_id=data.task_id,
                            # project_sn=data.project_sn,
                            project_sn=task_info['project_sn'],
                            module_type=module_type, params=params_json)
        post_info = super(AnnoHmdbAction, self).POST()
        if post_info['success']:
            post_info['content'] = {'ids': {'id': str(main_table_id), 'name': main_table_name}}
        return json.dumps(post_info)

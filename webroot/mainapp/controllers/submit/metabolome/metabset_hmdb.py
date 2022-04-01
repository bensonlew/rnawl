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
# __author__ = 'shaohua.yuan'


class MetabsetHmdbAction(MetabolomeController):
    def __init__(self):
        super(MetabsetHmdbAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        for arg in ['metabset', 'task_id', 'submit_location', 'task_type']:
            if not hasattr(data, arg):
                info = {'success': False, 'info': 'Argument:{} required'.format(arg)}
                return json.dumps(info)
        ## 和代谢集注释共用workflow
        task_name = 'metabolome.report.anno_hmdb'
        module_type = 'workflow'
        params_json = {
            'metabset': data.metabset,
            'task_id': data.task_id,
            'submit_location': data.submit_location,
            'task_type': int(data.task_type)
        }
        task_info = self.Metabolome.get_task_info(data.task_id)
        if task_info['type'] == 'GC':
            table_type = 'pos'
        else:
            table_type = 'mix'

        mongo_data = [
            # ('project_sn', data.project_sn),
            ('project_sn', task_info['project_sn']),
            ('task_id', data.task_id),
            ('status', 'start'),
            ('desc', '正在计算'),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ('metab_set',ObjectId(data.metabset)),
            ('table_type',table_type)
        ]
        metabolome = Metabolome()
        set_name_r = metabolome.conmon_find_one('metab_set',{"_id":ObjectId(data.metabset)})
        if set_name_r:
            set_name_pls = set_name_r['name']+'_'
        else:
            set_name_pls = ''
        main_table_name = "MetabsetHmdb_" +set_name_pls + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        mongo_data.append(('name', main_table_name))
        mongo_data.append(('params', json.dumps(params_json, sort_keys=True, separators=(',', ':'))))

        main_table_id = metabolome.insert_main_table('metabset_hmdb', mongo_data)
        metabolome.insert_main_id('metabset_hmdb', main_table_id)
        update_info = {str(main_table_id): 'metabset_hmdb'}
        hmdb_overview = metabolome.get_hmdb_anno(data.task_id)
        if hmdb_overview == "":
            info = {"success": False, "info": "this project has no origin hmdb annotation result" }
            return json.dumps(info)
        options = {
            'update_info': json.dumps(update_info),
            'anno_overview': hmdb_overview,
            'metabtable': data.metabset,
            'main_table_id': str(main_table_id),
            'name': main_table_name,
            "type": "metabsethmdb"
        }
        to_file = ['metabolome.export_metab_set(metabtable)']
        task_info = self.Metabolome.get_task_info(data.task_id)
        if "save_pdf" in task_info and int(task_info["save_pdf"]) == 1:
            options["save_pdf"] = 1
            m_table_name = "5.Metabset/08.MetabsetHMDB/" + main_table_name
        else:
            m_table_name = "Metabset/" + main_table_name.strip().split("_")[0] + '/' + main_table_name
        self.set_sheet_data(name=task_name, options=options,
                            main_table_name=m_table_name,
                            task_id=data.task_id,to_file=to_file,
                            project_sn=task_info['project_sn'],
                            module_type=module_type, params=params_json)
        post_info = super(MetabsetHmdbAction, self).POST()
        if post_info['success']:
            post_info['content'] = {'ids': {'id': str(main_table_id), 'name': main_table_name}}
        return json.dumps(post_info)

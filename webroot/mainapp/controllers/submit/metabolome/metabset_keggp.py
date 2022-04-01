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
# __author__ = 'guhaidong'


class MetabsetKeggpAction(MetabolomeController):
    def __init__(self):
        super(MetabsetKeggpAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        # for arg in ['metabset', 'task_id', 'submit_location','task_type','project_sn']:
        for arg in ['metabset', 'task_id', 'submit_location','task_type']:
            if not hasattr(data, arg):
                info = {'success': False, 'info': 'Argument:{} required'.format(arg)}
                return json.dumps(info)
        task_name = 'metabolome.report.metabset_keggp'
        module_type = 'workflow'
        params_json = {
            'metabset': data.metabset,
            'task_id': data.task_id,
            'task_type': int(data.task_type),
            'submit_location': data.submit_location,
            # 'project_sn': data.project_sn
        }
        task_info = self.Metabolome.get_task_info(data.task_id)
        mongo_data = [
            # ('project_sn', data.project_sn),
            ('project_sn', task_info['project_sn']),
            ('task_id', data.task_id),
            ('status', 'start'),
            ('desc', '正在计算'),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ('version','3.2')
        ]
        main_table_name = "MetabsetKeggp_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        mongo_data.append(('name', main_table_name))
        mongo_data.append(('params', json.dumps(params_json, sort_keys=True, separators=(',', ':'))))
        #rere_table_path = self._get_file_dir('/pathway_img', data.task_id, main_table_name.strip().split("_")[0] + '/' + main_table_name)
       # rere_table_path = rere_table_path + "/"
        rere_table_path = self._get_file_path('/pathway_img', data.task_id, "Metabset/"+main_table_name.strip().split("_")[0] + '/' + main_table_name)
        mongo_data.append(('result_dir', rere_table_path))
        metabolome = Metabolome()
        main_table_id = metabolome.insert_main_table('metabset_keggp', mongo_data)
        metabolome.insert_main_id('metabset_keggp', main_table_id)
        update_info = {str(main_table_id): 'metabset_keggp'}
        options = {
            'update_info': json.dumps(update_info),
            'metabset': data.metabset,
            'anno_overview': data.task_id,
            'trans_file': data.task_id,
            'main_table_id': str(main_table_id),
            'name': main_table_name
        }
        to_file = ['metabolome.export_mul_metab_set(metabset)', 'metabolome.export_overview_ko(anno_overview)', 'metabolome.export_overview(trans_file)']
        task_info = self.Metabolome.get_task_info(data.task_id)
        if "save_pdf" in task_info and int(task_info["save_pdf"]) == 1:
            options["save_pdf"] = 1
            m_table_name = "5.Metabset/06.MetabsetKeggp/" + main_table_name
        else:
            m_table_name = "Metabset/" + main_table_name.strip().split("_")[0] + '/' + main_table_name
        self.set_sheet_data(name=task_name, options=options,
                            main_table_name=m_table_name,
                            task_id=data.task_id,
                            # project_sn=data.project_sn,
                            project_sn=task_info['project_sn'],
                            module_type=module_type, params=params_json, to_file=to_file)
        post_info = super(MetabsetKeggpAction, self).POST()
        if post_info['success']:
            post_info['content'] = {'ids': {'id': str(main_table_id), 'name': main_table_name}}
        self.Metabolome.insert_set_info(data.metabset, "metabset_keggp", main_table_id)
        return json.dumps(post_info)

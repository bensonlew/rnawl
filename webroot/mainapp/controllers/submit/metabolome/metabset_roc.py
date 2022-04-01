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
# __author__ = 'zouguanqing'


class MetabsetRocAction(MetabolomeController):
    def __init__(self):
        super(MetabsetRocAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        print data
        # for arg in ['metabset', 'task_id', 'submit_location','task_type','project_sn']:
        for arg in ['metab_set', 'task_id', 'submit_location','task_type','group_id','group_detail','metab_table']:  #v3 rm 'compound'
            if not hasattr(data, arg):
                info = {'success': False, 'info': 'Argument:{} required'.format(arg)}
                return json.dumps(info)
        task_name = 'metabolome.report.metabset_roc'
        module_type = 'workflow'
        params_json = {
            'metab_set': data.metab_set,
            'task_id': data.task_id,
            'task_type': int(data.task_type),
            'submit_location': data.submit_location,
            'group_id' : data.group_id,
            'group_detail' : json.loads(data.group_detail),  #
            'metab_table' : data.metab_table
            #'compound' : data.compound
            # 'project_sn': data.project_sn
        }
        task_info = self.Metabolome.get_task_info(data.task_id)
        project_type = task_info['type']
        if project_type == 'LC':
            exp_table_path = self.Metabolome.get_metab_table(data.metab_table,data.task_id,type='mix')
            metab_desc_path = self.Metabolome.get_metab_desc(data.metab_table, data.task_id, type='mix')
        else:
            exp_table_path = self.Metabolome.get_metab_table(data.metab_table,data.task_id,type='pos')
            metab_desc_path = self.Metabolome.get_metab_desc(data.metab_table, data.task_id, type='pos')

        mongo_data = [
            # ('project_sn', data.project_sn),
            ('project_sn', task_info['project_sn']),
            ('task_id', data.task_id),
            ('status', 'start'),
            ('desc', '正在计算'),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ('version','3.0')
        ]
        main_table_name = "MetabsetRoc_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        mongo_data.append(('name', main_table_name))
        mongo_data.append(('params', json.dumps(params_json, sort_keys=True, separators=(',', ':'))))
        #rere_table_path = self._get_file_path('/pathway_img', data.task_id, main_table_name.strip().split("_")[0] + '/' + main_table_name)

        metabolome = Metabolome()
        main_table_id = metabolome.insert_main_table('metabset_roc', mongo_data)
        metabolome.insert_main_id('metabset_roc', main_table_id)
        update_info = {str(main_table_id): 'metabset_roc'}
        options = {
            'update_info': json.dumps(update_info),
            'exp_table' : exp_table_path,
            "metab_desc": metab_desc_path,
            'group_table' : data.group_id,
            "group_detail": data.group_detail,
            "metab_set_table": data.metab_set,
            #"compound" : data.compound,
            'roc_id': str(main_table_id),
            "name": main_table_name
        }

        to_file = []
        to_file.append('metabolome.export_metab_set1(metab_set_table)')
        to_file.append('metabolome.export_group_by_detail(group_table)')
        task_info = self.Metabolome.get_task_info(data.task_id)
        if "save_pdf" in task_info and int(task_info["save_pdf"]) == 1:
            options["save_pdf"] = 1
            m_table_name = "5.Metabset/10.MetabsetRoc/" + main_table_name
        else:
            m_table_name = "Metabset/" + main_table_name.strip().split("_")[0] + '/' + main_table_name
        self.set_sheet_data(name=task_name, options=options,
                            main_table_name=m_table_name,
                            task_id=data.task_id,
                            # project_sn=data.project_sn,
                            project_sn=task_info['project_sn'],
                            module_type=module_type, params=params_json, to_file=to_file)
        post_info = super(MetabsetRocAction, self).POST()
        if post_info['success']:
            post_info['content'] = {'ids': {'id': str(main_table_id), 'name': main_table_name}}
        self.Metabolome.insert_set_info(data.metab_set, "metabset_roc", main_table_id)
        return json.dumps(post_info)

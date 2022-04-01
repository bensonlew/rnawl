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
# __author__ = 'guhaidong'


class MetabsetIpathAction(MetabolomeController):
    def __init__(self):
        super(MetabsetIpathAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        # for arg in ['metabset', 'task_id', 'submit_location','task_type', 'project_sn']:
        for arg in ['metabset', 'task_id', 'submit_location','task_type']:
            if not hasattr(data, arg):
                info = {'success': False, 'info': 'Argument:{} required'.format(arg)}
                return json.dumps(info)
        version = 3
        if version ==2:
            task_name = 'metabolome.report.metabset_ipath'
        else:
            task_name = 'metabolome.report.metabset_ipath3'

        module_type = 'workflow'
        params_json = {
            'metabset': data.metabset,
            'task_id': data.task_id,
            'task_type': int(data.task_type),
            'submit_location': data.submit_location
            # 'project_sn': data.project_sn
        }
        task_info = self.Metabolome.get_task_info(data.task_id)
        mongo_data = [
            # ('project_sn', data.project_sn),
            ('project_sn', task_info['project_sn']),
            ('task_id', data.task_id),
            ('status', 'start'),
            ('desc', '正在计算'),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
             #####!!!
        ]
        if version ==2:
            mongo_data.append(("pathways", ["Metabolic_pathways.svg", "Regulatory_pathways.svg", "Biosynthesis_of_secondary_metabolities.svg"]))
        else:
            mongo_data.append(("pathways", ["Metabolism.svg", "Antibiotics.svg", "Microbial_metabolism.svg","Secondary_metabolites.svg"]))
        metabolome = Metabolome()
        # set_name_r = metabolome.conmon_find_one('metab_set',{"_id":ObjectId(data.metabset)})
        # if set_name_r:
        #     set_name_pls = set_name_r['name']+'_'
        # else:
        set_name_pls = ''
        main_table_name = "MetabsetIpath_"+set_name_pls + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        mongo_data.append(('name', main_table_name))
        mongo_data.append(('params', json.dumps(params_json, sort_keys=True, separators=(',', ':'))))
        #rere_table_path = self._get_file_dir('', data.task_id, main_table_name.strip().split("_")[0] + '/' + main_table_name)
        task_info = self.Metabolome.get_task_info(data.task_id)
        if "save_pdf" in task_info and int(task_info["save_pdf"]) == 1:
            # options["save_pdf"] = 1
            m_table_name = "5.Metabset/09.MetabsetIpath/" + main_table_name
        else:
            m_table_name = "Metabset/" + main_table_name.strip().split("_")[0] + '/' + main_table_name
        rere_table_path = self._get_file_path('',data.task_id, m_table_name)
        mongo_data.append(('result_dir', rere_table_path))

        main_table_id = metabolome.insert_main_table('metabset_ipath', mongo_data)
        metabolome.insert_main_id('metabset_ipath', main_table_id)
        update_info = {str(main_table_id): 'metabset_ipath'}
        options = {
            'metabset': data.metabset,
            'anno_overview': data.task_id,
            'main_table_id': str(main_table_id),
            'update_info': json.dumps(update_info),
            # "name": main_table_name
        }
        to_file = ['metabolome.export_mul_metab_set(metabset)', 'metabolome.export_overview(anno_overview)']
        self.set_sheet_data(name=task_name, options=options,
                            main_table_name=m_table_name,
                            task_id=data.task_id,
                            # project_sn=data.project_sn,
                            project_sn=task_info['project_sn'],
                            module_type=module_type, params=params_json, to_file=to_file)
        post_info = super(MetabsetIpathAction, self).POST()
        if post_info['success']:
            post_info['content'] = {'ids': {'id': str(main_table_id), 'name': main_table_name}}
        self.Metabolome.insert_set_info(data.metabset, "metabset_ipath", main_table_id)
        return json.dumps(post_info)

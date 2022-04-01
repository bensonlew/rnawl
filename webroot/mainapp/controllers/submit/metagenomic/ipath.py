# -*- coding: utf-8 -*-
# __author__ = 'guhaidong'

import web
import json
import datetime
from bson import ObjectId
from bson import SON
from mainapp.libs.param_pack import group_detail_sort
from mainapp.models.mongo.metagenomic import Metagenomic
from mainapp.controllers.project.metagenomic_controller import MetagenomicController
from mainapp.libs.signature import check_sig


class IpathAction(MetagenomicController):
    def __init__(self):
        super(IpathAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        default_argu = ['anno_id', 'geneset_id', 'group_detail', 'group_id',  'submit_location',  'task_type']
        for argu in default_argu:
            if not hasattr(data, argu):
                info = {"success": False, "info": "PARAMETERS MISSING: %s" , "variables": [argu], "code" : "C2403201"}
                return json.dumps(info)
        if not hasattr(data, 'anno_id'):
            info = {"success": False, "info": "PARAMETERS MISSING: anno_id", "code" : "C2403202"}
            return json.dumps(info)
        table_dict = json.loads(data.group_detail)
        if len(table_dict) < 2 or data.group_id == 'all':
            info = {"success": False, "info": "分析不适用于小于两个分组的情况！", "code":"C2403203", "variables":""}
            return  json.dumps(info)
        if not isinstance(table_dict, dict):
            info = {"success": False, "info": "PARAMETERS ERROR: wrong type of group_detail, dict expected!", "code" : "C2403204"}
            return  json.dumps(info)
        task_name = 'metagenomic.report.ipath'
        module_type = 'workflow'
        metagenomic = Metagenomic()
        geneset_info = metagenomic.get_geneset_info(data.geneset_id)
        task_info = metagenomic.get_task_info(geneset_info['task_id'])
        params_json = {
            'geneset_id': data.geneset_id,
            'group_id': data.group_id,
            'group_detail': group_detail_sort(data.group_detail),
            'task_type': int(data.task_type),
            'submit_location': data.submit_location,
            'anno_id': data.anno_id
        }
        #group_id = data.group_id if data.group_id in ['all', 'All', 'ALL'] else ObjectId(data.group_id)
        mongo_data = [
            ('project_sn', task_info['project_sn']),
            ('task_id', task_info['task_id']),
            ('status', 'start'),
            ('desc', 'processing'),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ("pathways", ["Metabolic_pathways.svg", "Regulatory_pathways.svg", "Biosynthesis_of_secondary_metabolities.svg"])
        ]
        main_table_name = 'Ipath_' + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        mongo_data.append(('name', main_table_name))
        mongo_data.append(('params', json.dumps(params_json, sort_keys=True, separators=(',', ':'))))
        main_table_id = self.metagenomic.insert_main_table('ipath', mongo_data)
        update_info = {str(main_table_id): 'ipath'}
        [geneset_table, gene_list] = metagenomic.export_geneset_table(data.geneset_id, "reads_num")
        options = {
            'update_info': json.dumps(update_info),
            'group_detail': data.group_detail,
            'group_file': data.geneset_id,
            'main_id': str(main_table_id),
            'geneset_table': geneset_table,
            "level_id": "KO",
            "color_style": "define1"
        }
        anno_info = metagenomic.get_anno_info(data.anno_id, "kegg")
        options['anno_table'] = anno_info['anno_file']
        options['lowest_level'] = anno_info['lowest_level']

        to_file = ['metagenomic.export_group_table_by_detail(group_file)']
        self.set_sheet_data(name=task_name, options=options,
                            main_table_name=main_table_name.strip().split("_")[0] + '/' + main_table_name,
                            task_id=task_info['task_id'],
                            project_sn=task_info['project_sn'], module_type=module_type, params=params_json,
                            to_file=to_file)
        task_info = super(IpathAction, self).POST()
        if task_info['success']:
            task_info['content'] = {'ids': {'id': str(main_table_id), 'name': main_table_name}}
        return json.dumps(task_info)

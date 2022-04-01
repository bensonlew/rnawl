# -*- coding: utf-8 -*-
# __author__ = 'gaohao'

import web
import json
import datetime
from bson import ObjectId
from bson import SON
from mainapp.libs.param_pack import group_detail_sort
from mainapp.models.mongo.metagenomic import Metagenomic
from mainapp.controllers.project.metagenomic_controller import MetagenomicController
from mainapp.libs.signature import check_sig


class MgRandomforestAction(MetagenomicController):
    def __init__(self):
        super(MgRandomforestAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        print data
        default_argu = ['anno_type', 'geneset_id', 'method', 'submit_location', 'group_detail','group_id', 'auth_method','tree_num']
        for argu in default_argu:
            if not hasattr(data, argu):
                info = {"success": False, "info": "PARAMETERS MISSING: %s" % argu, "code" : "C2402701"}
                return json.dumps(info)
        if data.auth_method not in [ 'CV','AUC']:
            info = {"success": False, "info": "PARAMETERS ERROR: wrong value of auth_method (%s)" % data.auth_method, "code" : "C2402702"}
            return json.dumps(info)
        if isinstance(data.tree_num,int):
            info = {"success": False, "info": "PARAMETERS ERROR: wrong value of tree_num (%s),is not int!" % data.tree_num, "code" : "C2402703"}
            return json.dumps(info)
        if data.anno_type not in ["kegg", "cog", "vfdb", "ardb", "card", "cazy","annopersonal","gene","nr"]:
            info = {"success": False, "info": "PARAMETERS ERROR: wrong value of anno_type (%s)" % data.anno_type, "code" : "C2402704"}
            return json.dumps(info)
        group_detail = json.loads(data.group_detail)
        if not isinstance(group_detail, dict):
            info = {"success": False, "info": 'PARAMETERS ERROR: wrong type of group_detail, dict expected!', "code" : "C2402705"}
            return json.dumps(info)
        if data.group_id == 'all':
            info = {"success": False, "info": 'Grouping scheme selects at least two groups！', "code" : "C2402706"}
            return json.dumps(info)
        elif len(group_detail) < 2:
            info = {"success": False, "info": 'Please select at least two groups of grouping schemes!', "code" : "C2402707"}
            return json.dumps(info)
        if hasattr(data, 'pre_group_id') and hasattr(data, 'pre_group_detail'):
            pre_group_detail = json.loads(data.pre_group_detail)
            if not isinstance(pre_group_detail, dict):
                info = {"success": False, "info": 'PARAMETERS ERROR: wrong type of pre_group_detail, dict expected!', "code" : "C2402708"}
                return json.dumps(info)
            if len(data.pre_group_detail) < 1:
                info = {"success": False, "info": 'PARAMETERS ERROR: wrong type of pre_group_detail!', "code" : "C2402709"}
                return json.dumps(info)
        else:
            pass
        key1 = list(group_detail.values())
        for i in range(len(key1)):
            if (len(key1[i]) < 5):
                info = {"success": False, "info": 'The number of samples in the group cannot be less than 5, please check!', "code" : "C2402710"}
                return json.dumps(info)
        if data.anno_type in ['nr']:
            if not hasattr(data, 'level_id'):
                info = {"success": False, "info": "PARAMETERS MISSING: level_id", "code" : "C2402711"}
                return json.dumps(info)
            if int(data.level_id) not in [1, 2, 3, 4, 5, 6, 7, 8, 9]:
                info = {"success": False, "info": 'PARAMETERS ERROR: wrong value of level_id', "code" : "C2402712"}
                return json.dumps(info)
        if data.anno_type in ['annopersonal']:
            if not hasattr(data, 'database'):
                info = {"success": False, "info": "PARAMETERS MISSING: database", "code" : "C2402713"}
                return json.dumps(info)
            if data.database not in ['go','phi','qs','mvirdb','tcdb','pfam','cyps','probio']:
                info = {"success": False, "info": "PARAMETERS ERROR: wrong value of database (%s)" % data.database, "code" : "C2402714"}
                return json.dumps(info)
        if data.anno_type not in ['nr','gene']:
            if not hasattr(data, 'level_id'):
                info = {"success": False, "info": "PARAMETERS MISSING: level_id", "code" : "C2402715"}
                return json.dumps(info)
        task_name = 'metagenomic.report.randomforest'
        module_type = 'workflow'
        task_type = 2
        metagenomic = Metagenomic()
        geneset_info = metagenomic.get_geneset_info(data.geneset_id)
        task_info = metagenomic.get_task_info(geneset_info['task_id'])
        group_id = data.group_id if data.group_id in ['all', 'All', 'ALL'] else ObjectId(data.group_id)
        params_json = {
            'anno_type': data.anno_type,
            'geneset_id': data.geneset_id,
            'method': data.method,
            'group_id': data.group_id,
            'group_detail': group_detail_sort(data.group_detail),
            'auth_method': data.auth_method,
            'task_type': int(task_type),
            'tree_num':int(data.tree_num),
            'submit_location': data.submit_location
        }
        mongo_data = [
            ('project_sn', task_info['project_sn']),
            ('task_id', task_info['task_id']),
            ('status', 'start'),
            ('desc', 'Randomforest 分析'),
            ('anno_type',data.anno_type),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        ]
        tax_level =''
        func_level=''
        if data.anno_type in ['nr']:
            params_json['anno_id'] = data.anno_id
            params_json['level_id'] = int(data.level_id)
            tax_level = self.level_id(data.level_id)
        if data.anno_type in ['annopersonal']:
            params_json['anno_id'] = data.anno_id
            params_json['database'] = data.database
            params_json['level_id'] = int(data.level_id)
            func_level = self.level_id(data.level_id)
        if data.anno_type not in ['gene','nr','annopersonal']:
            params_json['anno_id'] = data.anno_id
            params_json['level_id'] = int(data.level_id)
            func_level = self.level_id(data.level_id)
        if hasattr(data, 'pre_group_id') and hasattr(data, 'pre_group_detail'):
            params_json['pre_group_detail'] = group_detail_sort(data.pre_group_detail)
            params_json['pre_group_id'] = data.pre_group_id
        main_table_name = ''
        if data.anno_type in ['nr']:
            main_table_name = 'Randomforest' + '_' + data.anno_type.upper() + '_' + tax_level.capitalize() + '_' + \
                              datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        if data.anno_type in ['annopersonal']:
            main_table_name = 'Randomforest' +'_'+ data.database.upper() + '_'  + func_level.capitalize() + '_' + \
                          datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        if data.anno_type not in ['annopersonal','nr','gene']:
            main_table_name = 'Randomforest' + '_' + data.anno_type.upper() + '_' + func_level.capitalize() + '_' + \
                              datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        if data.anno_type in ['gene']:
            main_table_name = 'Randomforest' + '_' + data.anno_type.upper() + '_' +  \
                              datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        mongo_data.append(('name', main_table_name))
        mongo_data.append(('params', json.dumps(params_json, sort_keys=True, separators=(',', ':'))))
        main_table_id = self.metagenomic.insert_main_table('randomforest', mongo_data)
        update_info = {str(main_table_id): 'randomforest'}
        [geneset_table,gene_list]= metagenomic.export_geneset_table(data.geneset_id, data.method)
        to_file = []
        options = {
            'geneset_table': geneset_table,
            'update_info': json.dumps(update_info),
            'group_detail': data.group_detail,
            'group_id': data.group_id,
            'anno_type': data.anno_type,
            'auth_method': data.auth_method,
            'group': data.geneset_id,  # 为修改id转样品功能
            'tree_num': data.tree_num,
        }
        to_file.append('metagenomic.export_group_table_by_detail(group)')
        if data.anno_type in ['nr']:
            nr_anno_info = metagenomic.get_anno_info(data.anno_id, 'nr')
            nr_anno_table = nr_anno_info['anno_file']
            options['anno_table']=nr_anno_table
            options['level_type'] = self.level_id(data.level_id)
        if data.anno_type not in ['gene','nr']:
            if data.anno_type not in ['annopersonal']:
                func_anno_info = metagenomic.get_anno_info(data.anno_id, data.anno_type)
                func_anno_table = func_anno_info['anno_file']
            else:
                func_anno_info = metagenomic.get_anno_info(data.anno_id, data.database)
                func_anno_table = func_anno_info['anno_file']
            options['anno_table']=func_anno_table
            options['level_type'] = self.level_id(data.level_id)
            options['lowest_level'] = func_anno_info['lowest_level']
        if data.anno_type in ['gene']:
            options['gene_list'] = gene_list
        if hasattr(data, 'pre_group_id') and hasattr(data, 'pre_group_detail'):
            options['pre_group_detail'] = data.pre_group_detail
            options['pre_group_id'] = data.pre_group_id
            options['pre_group'] = data.geneset_id
            to_file.append('metagenomic.export_group_table_by_rf_detail(pre_group)')
        options['main_id'] = str(main_table_id)
        options['main_table_data'] = SON(mongo_data)
        self.set_sheet_data(name=task_name, options=options,
                            main_table_name=main_table_name.strip().split("_")[0] + '/' + main_table_name,
                            task_id=task_info['task_id'],
                            project_sn=task_info['project_sn'], module_type=module_type, params=params_json,
                            to_file=to_file)
        task_info = super(MgRandomforestAction, self).POST()
        if task_info['success']:
            task_info['content'] = {'ids': {'id': str(main_table_id), 'name': main_table_name}}
        return json.dumps(task_info)

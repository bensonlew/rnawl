# -*- coding: utf-8 -*-
# __author__ = 'zhujuan'

import web
import json
import datetime
from bson import ObjectId
from bson import SON
from mainapp.libs.param_pack import group_detail_sort
from mainapp.models.mongo.metagenomic import Metagenomic
from mainapp.controllers.project.metagenomic_controller import MetagenomicController
from mainapp.libs.signature import check_sig
from .comm_creat import CommCreatAction


class CompositionAction(MetagenomicController):
    def __init__(self):
        super(CompositionAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        print data
        default_argu = ['anno_type', 'geneset_id', 'method', 'submit_location', 'group_id', 'group_detail',
                        'group_method', 'graphic_type']
        for argu in default_argu:
            if not hasattr(data, argu):
                info = {"success": False, "info": "PARAMETERS MISSING: %s!" % argu, "code" : "C2401103"}
                return json.dumps(info)
        if data.group_method not in ["", "sum", "average", "middle"]:
            info = {"success": False, "info": "PARAMETERS ERROR: wrong value of group_method (%s)" % data.group_method, "code" : "C2401104"}
            return json.dumps(info)
        add_personal = CommCreatAction()
        info = add_personal.judge_database(data)
        if info != "":
            return json.dumps(info)
        task_name = 'metagenomic.report.composition'
        module_type = 'workflow'
        task_type = 2
        metagenomic = Metagenomic()
        geneset_info = metagenomic.get_geneset_info(data.geneset_id)
        task_info = metagenomic.get_task_info(geneset_info['task_id'])
        params_json = {
            'anno_type': data.anno_type,
            'geneset_id': data.geneset_id,
            'method': data.method,
            'group_id': data.group_id,
            'group_detail': group_detail_sort(data.group_detail),
            'group_method': data.group_method,
            'graphic_type': data.graphic_type,
            'task_type': task_type,
            'submit_location': data.submit_location
        }
        group_id = data.group_id if data.group_id in ['all', 'All', 'ALL'] else ObjectId(data.group_id)
        mongo_data = [
            ('project_sn', task_info['project_sn']),
            ('task_id', task_info['task_id']),
            #('geneset_id', ObjectId(data.geneset_id)),
            ('status', 'start'),
            #('group_id', group_id),
            ('desc', 'processing'),
            ('anno_type', data.anno_type),
            ('graphic_type', data.graphic_type),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        ]
        level = ""
        level_name = ""
        if data.anno_type not in ["gene"]:
            params_json['anno_id'] = data.anno_id
            params_json['level_id'] = int(data.level_id)  # 此处是level缩写，需要再转义成真实level name
            level = data.level_id
            level_name = self.level_id(data.level_id)
            #mongo_data.append(('anno_id', ObjectId(data.anno_id)))
            #mongo_data.append(('level_id', data.level_id))
            if hasattr(data, 'second_level'):
                params_json['second_level'] = data.second_level
                # mongo_data.append(('second_level', data.second_level))
        if data.graphic_type in ["bar", "circos"]:
            if hasattr(data, 'combine_value'):
                params_json['combine_value'] = float(data.combine_value) if float(data.combine_value) !=0 else 0
            else:
                info = {'success': False, 'info': '参数缺失：Others合并', 'code':'C2401101', 'variables':''}
                return json.dumps(info)
        elif data.graphic_type in ["heatmap"]:
            if hasattr(data, 'top'):
                params_json['top'] = int(data.top)
            else:
                info = {'success': False, 'info': '参数缺失：分类水平总丰度前N', 'code':'C2401102', 'variables':''}
                return json.dumps(info)
            if hasattr(data, 'species_cluster'):
                params_json['species_cluster'] = data.species_cluster
            else:
                info = {'success': False, 'info': 'PARAMETERS MISSING: species_cluster', "code" : "C2401108"}
                return json.dumps(info)
            if hasattr(data, 'specimen_cluster'):
                params_json['specimen_cluster'] = data.specimen_cluster
            else:
                info = {'success': False, 'info': 'PARAMETERS MISSING: specimen_cluster', "code" : "C2401109"}
                return json.dumps(info)
            if hasattr(data, "color_level"):
                params_json["color_level"]= int(data.color_level)#add by qingchen.zhang@20181202用于颜色水平的筛选
            if hasattr(data, "normalization"):
                params_json["normalization"] = data.normalization
        elif data.graphic_type in ["bubble"]:
            if hasattr(data, 'combine_value'):
                params_json['combine_value'] = float(data.combine_value)
            else:
                info = {'success': False, 'info': '参数缺失：Others合并', 'code':'C2401101', 'variables':''}
                return json.dumps(info)
            if hasattr(data, "color_level"):
                params_json["color_level"]= int(data.color_level) #add by qingchen.zhang@20181202用于颜色水平的筛选
        [geneset_table, gene_list] = metagenomic.export_geneset_table(data.geneset_id, data.method)
        to_file = []
        options = {
            'graphic_type': data.graphic_type,
            'method': data.method,
            'group_method': data.group_method,
            #'update_info': json.dumps(update_info),
            'params': json.dumps(params_json, sort_keys=True, separators=(',', ':')),
            'group_id': data.group_id,
            'geneset_id': data.geneset_id,
            'group_detail': data.group_detail,
            'geneset_table': geneset_table,
            'anno_type': data.anno_type,
            'group': data.geneset_id,  # 为修改id转样品功能
            #'project_name': 'metagenomic',#为与meta的项目区分开
        }
        #个性化数据库调用
        add_personal = CommCreatAction()
        data,options, params_json,info, to_file = add_personal.add_personal(data,options,params_json,geneset_table,to_file)
        if info!= "":
            return json.dumps(info)
        #个性化数据库用调用
        to_file.append('metagenomic.export_group_table_by_detail(group)')
        if data.anno_type not in ['gene']:
            options['level_id'] = self.level_id(data.level_id)
            anno_info = metagenomic.get_anno_info(data.anno_id, data.anno_type)
            options['anno_id'] = data.anno_id
            options['anno_table'] = anno_info['anno_file']
            if data.anno_type != 'nr':
                options['lowest_level'] = anno_info['lowest_level']
            if hasattr(data, 'second_level'):
                if data.second_level != "":
                    options['level_type_name'] = self.level_convert(data.second_level, int(data.level_id))
        else:
            options['gene_list'] = gene_list
        if hasattr(data, "normalization"):
            options["normalization"] = '' if data.normalization == "none" else data.normalization
        if hasattr(data, 'combine_value'):
            options['others'] = float(data.combine_value)
        if hasattr(data, 'top'):
            options['species_number'] = data.top
        if hasattr(data, 'specimen_cluster'):
            options['sample_method'] = data.specimen_cluster
        if hasattr(data, 'species_cluster'):
            options['species_method'] = data.species_cluster
        if hasattr(data, "color_level"):
            options['level_color'] = self.level_id(data.color_level)
        main_table_name = data.graphic_type.capitalize() + '_' + data.anno_type.upper() + '_' + level_name + '_' + \
                          datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        main_table_name = main_table_name.replace(" ","_")
        mongo_data.append(('name', main_table_name))
        mongo_data.append(('params', json.dumps(params_json, sort_keys=True, separators=(',', ':'))))
        main_table_id = self.metagenomic.insert_main_table('composition', mongo_data)
        update_info = {str(main_table_id): 'composition'}
        options['main_id'] = str(main_table_id)
        options['main_table_data'] = SON(mongo_data)
        options['update_info'] = json.dumps(update_info)
        self.set_sheet_data(name=task_name, options=options,
                            main_table_name=main_table_name.strip().split("_")[0] + '/' + main_table_name,
                            task_id=task_info['task_id'],
                            project_sn=task_info['project_sn'], module_type=module_type, params=params_json,
                            to_file=to_file)
        task_info = super(CompositionAction, self).POST()
        if task_info['success']:
            task_info['content'] = {'ids': {'id': str(main_table_id), 'name': main_table_name}}
        return json.dumps(task_info)

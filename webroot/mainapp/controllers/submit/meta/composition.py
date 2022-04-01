# -*- coding: utf-8 -*-
# __author__ = 'zhujuan'


import web
import json
import datetime
from bson import ObjectId
from bson import SON
from mainapp.libs.param_pack import group_detail_sort,param_pack
from mainapp.libs.signature import check_sig
from mainapp.controllers.project.meta_controller import MetaController

class CompositionAction(MetaController):
    def __init__(self):
        super(CompositionAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        #default_argu = ['anno_type', 'geneset_id', 'method', 'submit_location', 'group_id', 'group_detail',
        #                'group_method', 'graphic_type']
        default_argu = ['otu_id', 'submit_location', 'group_id', 'group_detail', 'group_method', 'graphic_type']
        for argu in default_argu:
            if not hasattr(data, argu):
                info = {"success": False, "info": "parameters missing:%s" % argu}
                return json.dumps(info)
        if data.group_method not in ["", "sum", "average", "middle"]:
            variables = []
            variables.append(data.group_method)
            info = {"success": False, "info": "对分组样本计算方式:%s错误!" % data.group_method, 'code':'C2200501', 'variables':variables}
            return json.dumps(info)
        #if data.anno_type not in ["nr", "kegg", "cog", "vfdb", "ardb", "card", "cazy", "gene"]:
        #    info = {"success": False, "info": "数据库类型:%s错误!" % data.anno_type}
        #    return json.dumps(info)
        #if data.anno_type not in ["gene"]:
        #    if not hasattr(data, 'anno_id'):
        #        info = {"success": False, "info": "缺少anno_id参数!"}
        #        return json.dumps(info)
        if not hasattr(data, 'level_id'):
            info = {"success": False, "info": "缺少level_id参数!", 'code':'C2200502'}
            return json.dumps(info)

        task_name = 'meta.report.composition'
        module_type = 'workflow'
        task_type = 1
        #metagenomic = Metagenomic()
        #geneset_info = metagenomic.get_geneset_info(data.geneset_id)
        otu_info = self.meta.get_otu_table_info(data.otu_id)
        if not otu_info:
            info = {"success": False, "info": "OTU不存在，请确认参数是否正确！", 'code':'C2200503'}
            return json.dumps(info)
        task_info = self.meta.get_task_info(otu_info['task_id'])

        #task_info = metagenomic.get_task_info(geneset_info['task_id'])
        submit_location = data.submit_location
        if data.graphic_type in ['ternary','Ternary']:  #guanqing.zou 20180511  页面记录区分circos
            submit_location = 'ternary'
        params_json = {
            'group_id': data.group_id,
            'group_detail': group_detail_sort(data.group_detail),
            'group_method': data.group_method,
            'graphic_type': data.graphic_type,
            'task_type': "reportTask", #guanqing.zou 20180503
            'submit_location' : submit_location,
            'combine_value' : data.combine_value,
            'level_id' : data.level_id,     #guanqing.zou 20180503
            'otu_id':data.otu_id              #guanqing.zou 20180503


        }
        group_id = data.group_id if data.group_id in ['all', 'All', 'ALL'] else ObjectId(data.group_id)
        mongo_data = [
            ('otu_id',ObjectId(data.otu_id)),
            ('from_id',data.otu_id),
            ('newick_id',None),
            ('show',0),
            ('project_sn', task_info['project_sn']),
            ('task_id', task_info['task_id']),
            #('geneset_id', ObjectId(data.geneset_id)),
            ('status', 'start'),
            #('group_id', group_id),
            ('desc', '正在计算'),
            #('anno_type', data.anno_type),
            ('type', data.graphic_type),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),

        ]
        level = ""
        level_name = ""
        # if data.anno_type not in ["gene"]:
        #    params_json['anno_id'] = data.anno_id
        #   params_json['level_id'] = int(data.level_id)  # 此处是level缩写，需要再转义成真实level name
        #   level = data.level_id
        #    level_name = self.level_id(data.level_id)
        #    #mongo_data.append(('anno_id', ObjectId(data.anno_id)))
        #    #mongo_data.append(('level_id', data.level_id))
        #    if hasattr(data, 'second_level'):
        #        params_json['second_level'] = data.second_level
        #       # mongo_data.append(('second_level', data.second_level))
        main_table_name = data.graphic_type.capitalize() + '_' + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        mongo_data.append(('name', main_table_name))
        if data.graphic_type in ["ternary", "Ternary"]:
            if hasattr(data, 'color_level'):
                params_json['color_level'] = data.color_level
                mongo_data.append(('color_level',data.color_level))
            else:
                info = {'success': False, 'info': '没有提供color_level', 'code':'C2200504'}
                return json.dumps(info)
        #if data.graphic_type in ["heatmap"]:
        #    if hasattr(data, 'top'):
        #        params_json['top'] = int(data.top)
        #    else:
        #        info = {'success': False, 'info': '没有提供选择分类水平总丰度前'}
        #        return json.dumps(info)
        #    if hasattr(data, 'species_cluster'):
        #        params_json['species_cluster'] = data.species_cluster
        #    else:
        #        info = {'success': False, 'info': '没有提供物种聚类方式'}
        #        return json.dumps(info)
        #    if hasattr(data, 'specimen_cluster'):
        #        params_json['specimen_cluster'] = data.specimen_cluster
        #    else:
        #        info = {'success': False, 'info': '没有提供样品聚类方式'}
        #        return json.dumps(info)
        mongo_data.append(('params', json.dumps(params_json, sort_keys=True, separators=(',', ':'))))

        #main_table_id = self.metagenomic.insert_main_table('composition', mongo_data)
        main_table_id = self.meta.insert_none_table('sg_composition')
        #update_info = {str(main_table_id): 'composition'}
        update_info = {str(main_table_id): 'sg_composition'}
        #[geneset_table, gene_list] = metagenomic.export_geneset_table(data.geneset_id, data.method)
        #to_file = []
        #options = {
        #    'type': data.graphic_type,
        #    'group_method': data.group_method,
        #    'update_info': json.dumps(update_info),
        #    'params': json.dumps(params_json, sort_keys=True, separators=(',', ':')),
        #    'group_id': data.group_id
            #'group': data.geneset_id  # 为修改id转样品功能
        #}
        
        options = {
            "input_otu_id": data.otu_id,
            "in_otu_table": data.otu_id,
            "group_detail": data.group_detail,
            "level": str(data.level_id),
            'update_info': json.dumps(update_info),
            'main_id': str(main_table_id),
            'main_table_data': SON(mongo_data),
            'method' : data.group_method,
            'combine_value' : data.combine_value,
            'type' : data.graphic_type
            }
        #to_file.append('metagenomic.export_group_table_by_detail(group)')
        to_file = "meta.export_otu_table_by_level(in_otu_table)"

        '''
        if data.anno_type not in ['gene']:
            options['level_id'] = self.level_id(data.level_id)
            anno_info = metagenomic.get_anno_info(data.anno_id, data.anno_type)
            options['anno_id'] = data.anno_id
            options['anno_table'] = anno_info['anno_file']
            if data.anno_type not in 'nr':
                options['lowest_level'] = anno_info['lowest_level']
            if hasattr(data, 'second_level'):
                if data.second_level != "":
                    options['level_type_name'] = self.level_convert(data.second_level, data.level_id)
        else:
            options['gene_list'] = gene_list
        if hasattr(data, 'combine_value'):
            options['others'] = data.combine_value
        if hasattr(data, 'top'):
            options['species_number'] = data.top
        if hasattr(data, 'specimen_cluster'):
            options['sample_method'] = data.specimen_cluster
        if hasattr(data, 'species_cluster'):
            options['species_method'] = data.species_cluster
        '''
        #options['main_id'] = str(main_table_id)
        #options['main_table_data'] = SON(mongo_data)
        self.set_sheet_data(name=task_name, options=options,
                            main_table_name=data.graphic_type.capitalize() + '/' + main_table_name,
                #           task_id=task_info['task_id'],project_sn=task_info['project_sn'], 
                            params=params_json,
                            module_type='workflow', to_file=to_file)
        task_info = super(CompositionAction, self).POST()
        if task_info['success']:
            task_info['content'] = {'ids': {'id': str(main_table_id), 'name': main_table_name}}
        return json.dumps(task_info)

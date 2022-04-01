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

class HeatmapCorAction(MetagenomicController):
    def __init__(self):
        super(HeatmapCorAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        default_argu = ['anno_type', 'geneset_id', 'method', 'submit_location', 'group_id',  'group_detail',
                        'env_id', 'env_labs', 'top', 'coefficient', 'species_cluster', 'env_cluster']
        for argu in default_argu:
            if not hasattr(data, argu):
                info = {"success": False, "info": "PARAMETERS MISSING: %s" % argu}
                return json.dumps(info)
        add_personal = CommCreatAction()
        info = add_personal.judge_database(data)
        if info != "":
            return json.dumps(info)
        if len(str(data.env_labs).split(",")) < 2:
            info = {'success': False, 'info': '相关性Heatmap分析环境因子数量必须大于等于2!', 'code':'C2401601', 'variables':''}
            return json.dumps(info)
        task_name = 'metagenomic.report.heatmap_cor'
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
            'submit_location': data.submit_location,
            'task_type': task_type,
            'env_id': data.env_id,
            'env_labs': data.env_labs,
            'env_cluster': data.env_cluster,
            'species_cluster': data.species_cluster,
            'top': int(data.top),
            'coefficient': data.coefficient
        }
        #group_id = data.group_id if data.group_id in ['all', 'All', 'ALL'] else ObjectId(data.group_id)
        mongo_data = [
            ('project_sn', task_info['project_sn']),
            ('task_id', task_info['task_id']),
            #('geneset_id', ObjectId(data.geneset_id)),
            ('status', 'start'),
            #('group_id', group_id),
            ('desc', 'processing'),
            ('anno_type', data.anno_type),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            #("env_id", ObjectId(data.env_id)),
            ('env_list', data.env_labs)
        ]
        #level = ""
        level_name = ""
        if data.anno_type not in ["gene"]:
            params_json['anno_id'] = data.anno_id
            params_json['level_id'] = int(data.level_id)
            #level = data.level_id
            level_name = self.level_id(data.level_id)
            #mongo_data.append(('anno_id', ObjectId(data.anno_id)))
            #mongo_data.append(('level_id', data.level_id))
            if hasattr(data, 'second_level'):
                params_json['second_level'] = data.second_level
                #mongo_data.append(('second_level', data.second_level))
        if data.coefficient == "pearsonr":
            method_name = "Pearson"
        elif data.coefficient == "spearmanr":
            method_name = "Spearman"
        else:
            method_name = "Kendall"
        #mongo_data.append(('name', main_table_name))
        #mongo_data.append(('params', json.dumps(params_json, sort_keys=True, separators=(',', ':'))))
        #main_table_id = self.metagenomic.insert_main_table('heatmap_cor', mongo_data)
        #update_info = {str(main_table_id): 'heatmap_cor'}
        [geneset_table, gene_list] = metagenomic.export_geneset_table(data.geneset_id, data.method)
        to_file = []
        options = {
            'env_id': data.env_id,
            "env_file": data.env_id,
            #'update_info': json.dumps(update_info),
            'env_labs': data.env_labs,
            "main_table_data": SON(mongo_data),
            'params': json.dumps(params_json, sort_keys=True, separators=(',', ':')),
            'group_id': data.group_id,
            'geneset_id': data.geneset_id,
            'group_detail': data.group_detail,
            'geneset_table': geneset_table,
            'anno_type': data.anno_type,
            'group': data.geneset_id,  # 为修改id转样品功能
            'top_species': data.top,
            'method': data.method,
            'analysis_method': data.coefficient,
            'env_cluster': data.env_cluster,
            'species_cluster': data.species_cluster,
        }
        #个性化数据库调用
        add_personal = CommCreatAction()
        data,options, params_json,info, to_file = add_personal.add_personal(data,options,params_json,geneset_table,to_file)
        if info!= "":
            return json.dumps(info)
        #个性化数据库用调用
        to_file.append('metagenomic.export_group_table_by_detail(group)')
        to_file.append('metagenomic.export_float_env(env_file)')
        if data.anno_type not in ['gene']:
            options['level_id'] = self.level_id(data.level_id)
            anno_info = metagenomic.get_anno_info(data.anno_id, data.anno_type)
            options['anno_id'] = data.anno_id
            options['anno_table'] = self.use_s3(anno_info['anno_file'])
            if data.anno_type not in 'nr':
                options['lowest_level'] = self.use_s3(anno_info['lowest_level'])
            if hasattr(data, 'second_level'):
                if data.second_level != "":
                    options['second_level'] = self.level_convert(data.second_level, data.level_id)
        else:
            options['gene_list'] = gene_list
        #    to_file.append('metagenomic.export_geneset_list(gene_list)')
        main_table_name = 'CorrHeatmap' + method_name + '_' + data.anno_type.upper() + '_' + level_name + '_' + \
            datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        main_table_name = main_table_name.replace(" ","_")
        mongo_data.append(('name', main_table_name))
        mongo_data.append(('params', json.dumps(params_json, sort_keys=True, separators=(',', ':'))))
        main_table_id = self.metagenomic.insert_main_table('heatmap_cor', mongo_data)
        update_info = {str(main_table_id): 'heatmap_cor'}
        options['main_id'] = str(main_table_id)
        options['main_table_data'] = SON(mongo_data)
        options['update_info'] = json.dumps(update_info)
        self.set_sheet_data(name=task_name, options=options,
                            main_table_name=main_table_name.strip().split("_")[0] + '/' + main_table_name,
                            task_id=task_info['task_id'], project_sn=task_info['project_sn'],
                            module_type=module_type, params=params_json, to_file=to_file)
        task_info = super(HeatmapCorAction, self).POST()
        if task_info['success']:
            task_info['content'] = {'ids': {'id': str(main_table_id), 'name': main_table_name}}
        return json.dumps(task_info)

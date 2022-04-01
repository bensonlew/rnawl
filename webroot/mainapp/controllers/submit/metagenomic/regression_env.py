
# -*- coding: utf-8 -*-
# __author__ = 'zouguanqing'

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


class RegressionEnvAction(MetagenomicController):
    def __init__(self):
        super(RegressionEnvAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        default_argu = ['anno_type', 'geneset_id', 'method', 'submit_location', 'group_detail','group_id', 'diversity_analysis_type','env_id','env_labs']
        for argu in default_argu:
            if not hasattr(data, argu):
                info = {"success": False, "info": "PARAMETERS MISSING: %s" , "variables": argu, "code" : "C2404401"}
                return json.dumps(info)
        if len(data.env_labs.split(',')) > 1:
            info = {"success": False, "info": "Select ONE ENV", "code" : "C2404402"}
            return json.dumps(info)

        if data.diversity_analysis_type not in ['pca','pcoa','nmds']:
            info = {"success": False, "info": "PARAMETERS ERROR: wrong value of diversity_analysis_type (%s)" , "variables": data.diversity_analysis_type, "code" : "C2404403"}
            return json.dumps(info)
        if data.anno_type not in ["kegg", "cog", "vfdb", "ardb", "card", "cazy", "nr", "gene", "annopersonal"]:    # add nr
            info = {"success": False, "info": "PARAMETERS ERROR: wrong value of anno_type (%s)" , "variables": data.anno_type, "code" : "C2404404"}
            return json.dumps(info)

        if data.anno_type != 'gene':
            if not hasattr(data, 'func_anno_id'):
                info = {"success": False, "info": "PARAMETERS MISSING: func_anno_id", "code" : "C2404405"}
                return json.dumps(info)
            if not hasattr(data,'tax_level_id'):
                info = {"success": False, "info": "PARAMETERS MISSING: tax_level_id", "code" : "C2404406"}
                return json.dumps(info)

        if data.diversity_analysis_type in ['pcoa', 'nmds']:
            if not hasattr(data, 'distance_type'):
                info = {"success": False, "info": "PARAMETERS MISSING: distance_type", "code" : "C2404407"}
                return json.dumps(info)
        task_name = 'metagenomic.report.regression_env'
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
            #'diversity_type': data.diversity_type,
            'diversity_analysis_type': data.diversity_analysis_type,
            #'tax_level_id': int(data.tax_level_id),
            'task_type': task_type,
            #'nr_anno_id': data.nr_anno_id,
            'submit_location': data.submit_location,
            'env_id': data.env_id,
            'env_labs': data.env_labs,
        }
        if data.anno_type != 'gene':
            params_json['tax_level_id'] = int(data.tax_level_id)
            params_json['func_anno_id'] = data.func_anno_id

        #group_id = data.group_id if data.group_id in ['all', 'All', 'ALL'] else ObjectId(data.group_id)
        if data.diversity_analysis_type != 'nmds':
            lab = 'PC1'
        else:
            lab = 'MDS1'
        mongo_data = [
            ('project_sn', task_info['project_sn']),
            ('task_id', task_info['task_id']),
            ('status', 'start'),
            ('desc', 'processing'),
            ('anno_type', data.anno_type),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ('table_type','env'),
            ('pcs', lab)
        ]
        #tax_level = self.level_id(data.tax_level_id

        if data.anno_type != 'gene':
            func_level = self.level_id(data.tax_level_id)
            #main_table_name = 'RegressionEnv_' + data.anno_type.upper() + '_'  + func_level.capitalize() + '_' + \
            #              datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]   #### 加 env id
        else:
            main_table_name = 'RegressionEnv_' + data.anno_type.upper() + '_' + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]   #### 加 env id

        if data.diversity_analysis_type in ['pcoa','nmds']:
            params_json['distance_type'] = data.distance_type
        [geneset_table,gene_list]= metagenomic.export_geneset_table(data.geneset_id, data.method)

        to_file = []
        options = {
            'geneset_table':geneset_table,
            #'update_info': json.dumps(update_info),
            'group_detail': data.group_detail,
            'group_id': data.group_id,
            'geneset_id': data.geneset_id,
            'anno_type': data.anno_type,
            'diversity_analysis_type': data.diversity_analysis_type,
            'group': data.geneset_id, # 为修改id转样品功能
            'envtable' : data.env_id,
            'env_labs' : data.env_labs
        }
        to_file.append('metagenomic.export_group_table_by_detail(group)')
        to_file.append('metagenomic.export_env_table_select(envtable)')
        #个性化数据库调用
        add_personal = CommCreatAction()
        data,options, params_json,info, to_file = add_personal.add_personal(data,options,params_json,geneset_table,to_file)
        if info!= "":
            return json.dumps(info)
        #个性化数据库用调用
        if data.anno_type not in ['gene']:
            options['level_id'] = func_level
            func_anno_info = metagenomic.get_anno_info(data.func_anno_id, data.anno_type)
            func_anno_table = func_anno_info['anno_file']
            options['anno_table'] = func_anno_table
            if data.anno_type not in 'nr':
                options['lowest_level'] = func_anno_info['lowest_level']
            if hasattr(data, 'second_level'):
                params_json['second_level'] =  data.second_level
                if data.second_level != "":
                    options['level_type_name'] = self.level_convert(data.second_level, data.tax_level_id)
            main_table_name = 'RegressionEnv_' + data.anno_type.upper() + '_'  + func_level.capitalize() + '_' + \
                          datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]   #### 加 env id
            main_table_name = main_table_name.replace(" ","_")
        else:
            options['gene_list'] = gene_list

        if data.diversity_analysis_type in ['pcoa', 'nmds']:
            options['distance_type'] = data.distance_type
        mongo_data.append(('name', main_table_name))
        mongo_data.append(('params', json.dumps(params_json, sort_keys=True, separators=(',', ':'))))
        main_table_id = self.metagenomic.insert_main_table('regression', mongo_data)
        update_info = {str(main_table_id): 'regression'}
        options['main_id'] = str(main_table_id)
        options['main_table_data'] = SON(mongo_data)
        options['update_info'] = json.dumps(update_info)
        self.set_sheet_data(name=task_name, options=options,
                            main_table_name=main_table_name.strip().split("_")[0] + '/' + main_table_name,
                            task_id=task_info['task_id'],
                            project_sn=task_info['project_sn'], module_type=module_type, params=params_json,
                            to_file=to_file)
        task_info = super(RegressionEnvAction, self).POST()
        if task_info['success']:
            task_info['content'] = {'ids': {'id': str(main_table_id), 'name': main_table_name}}
        return json.dumps(task_info)

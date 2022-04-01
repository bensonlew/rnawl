
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


class VpaAction(MetagenomicController):
    def __init__(self):
        super(VpaAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        print(data)
        default_argu = ['anno_type', 'geneset_id', 'method', 'submit_location', 'group_detail','task_type',
                         'group_id', 'env_id','env_group_id', 'env_detail']
        for argu in default_argu:
            if not hasattr(data, argu):
                info = {"success": False, "info": "PARAMETERS MISSING: %s" % argu}
                return json.dumps(info)
        add_personal = CommCreatAction()
        info = add_personal.judge_database(data)
        if info != "":
            return json.dumps(info)
        task_name = 'metagenomic.report.vpa'
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
            'task_type': task_type,
            'submit_location': data.submit_location,
            'env_id': data.env_id,
            #'env_labs': data.env_labs
            'env_detail': group_detail_sort(data.env_detail),
            'env_group_id': data.env_group_id
        }

        func_level = ''
        if data.anno_type not in ["gene"]:
            params_json['anno_id'] = data.anno_id
            params_json['level_id'] = int(data.level_id)
            func_level = self.level_id(data.level_id)
            if hasattr(data, 'second_level'):
                params_json['second_level'] = data.second_level


        #group_id = data.group_id if data.group_id in ['all', 'All', 'ALL'] else ObjectId(data.group_id)
        mongo_data = [
            ('project_sn', task_info['project_sn']),
            ('task_id', task_info['task_id']),
            ('status', 'start'),
            ('desc', 'processing'),
            ('anno_type', data.anno_type),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        ]

        [geneset_table,gene_list]= metagenomic.export_geneset_table(data.geneset_id, data.method)
        to_file = []
        options = {
            'geneset_table':geneset_table,
            #'update_info': json.dumps(update_info),
            'group_detail': data.group_detail,
            'params': json.dumps(params_json, sort_keys=True, separators=(',', ':')),
            'group_id': data.group_id,
            'geneset_id': data.geneset_id,
            'anno_type': data.anno_type,
            'env_detail' : data.env_detail,
            'group': data.geneset_id,
            'env_file': data.env_id,
            'env_id': data.env_id,
            'env_group':data.env_detail,
            #'main_id' : str(main_table_id)
        }
        #个性化数据库调用
        add_personal = CommCreatAction()
        data,options, params_json,info, to_file = add_personal.add_personal(data,options,params_json,geneset_table,to_file)
        if info!= "":
            return json.dumps(info)
        #个性化数据库用调用
        if data.anno_type not in ['gene']:
            options['level_id'] = func_level
            func_anno_info = metagenomic.get_anno_info(data.anno_id, data.anno_type)
            func_anno_table = func_anno_info['anno_file']
            options['anno_table'] = func_anno_table
            if data.anno_type not in 'nr':
                options['lowest_level'] = func_anno_info['lowest_level']
            if hasattr(data, 'second_level'):
                if data.second_level != "":
                    options['level_type_name'] = self.level_convert(data.second_level, data.level_id)
        else:
            options['gene_list'] = gene_list
        main_table_name = 'VPA_' + data.anno_type.upper() + '_'  + func_level.capitalize() + '_' + \
                          datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]   #### 加 env id
        main_table_name = main_table_name.replace(" ","_")
        mongo_data.append(('name', main_table_name))
        mongo_data.append(('params', json.dumps(params_json, sort_keys=True, separators=(',', ':'))))
        main_table_id = self.metagenomic.insert_main_table('vpa', mongo_data)
        update_info = {str(main_table_id): 'vpa'}
        to_file.append('metagenomic.export_group_table_by_detail(group)')
        to_file.append('metagenomic.export_env_table(env_file)')
        to_file.append('metagenomic.export_env_table_by_detail(env_group)')
        options['main_table_data'] = SON(mongo_data)
        options['update_info'] = json.dumps(update_info)
        options['main_id'] = str(main_table_id)
        self.set_sheet_data(name=task_name, options=options,
                            main_table_name=main_table_name.strip().split("_")[0] + '/' + main_table_name,
                            task_id=task_info['task_id'],
                            project_sn=task_info['project_sn'], module_type=module_type, params=params_json,
                            to_file=to_file)
        task_info = super(VpaAction, self).POST()
        if task_info['success']:
            task_info['content'] = {'ids': {'id': str(main_table_id), 'name': main_table_name}}
        return json.dumps(task_info)

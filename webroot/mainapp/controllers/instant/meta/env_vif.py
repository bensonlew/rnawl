# -*- coding: utf-8 -*-
# __author__ = 'zhujuan'

import web
import json
import datetime
from bson import ObjectId
from bson import SON
from mainapp.libs.param_pack import group_detail_sort
from mainapp.libs.signature import check_sig
from mainapp.controllers.project.meta_controller import MetaController


class EnvVifAction(MetaController):
    def __init__(self):
        super(EnvVifAction, self).__init__(instant=True)

    @check_sig
    def POST(self):
        data = web.input()
        default_argu = ['otu_id','level_id','submit_location', 'group_id', 'group_detail',"env_id", "env_labs", "viflim"]
        for argu in default_argu:
            if not hasattr(data, argu):
                info = {"success": False, "info": "parameters missing:%s" % argu}
                return json.dumps(info)

        task_name = 'meta.report.env_vif'
        module_type = 'workflow'
        task_type = 1

        otu_info = self.meta.get_otu_table_info(data.otu_id)
        task_info = self.meta.get_task_info(otu_info['task_id'])

        params_json = {
            'otu_id': data.otu_id,
            'level_id': int(data.level_id),
            'group_id': data.group_id,
            'group_detail': group_detail_sort(data.group_detail),
            'env_id': data.env_id,
            'env_labs': data.env_labs,
            'task_type': task_type,
            'submit_location': data.submit_location,
            'viflim': int(data.viflim)
        }
        group_id = data.group_id if data.group_id in ['all', 'All', 'ALL'] else ObjectId(data.group_id)
        mongo_data = [
            ('project_sn', task_info['project_sn']),
            ('task_id', task_info['task_id']),
            ('otu_id',ObjectId(data.otu_id)),  
            ('from_id',data.otu_id),
            ('status', 'start'),
            # ('group_id', group_id),
            # ('env_id', ObjectId(data.env_id)),
            ('desc', '正在计算'),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        ]

        main_table_name = "Vif" + '_'  +  datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        mongo_data.append(('name', main_table_name))
        mongo_data.append(('params', json.dumps(params_json, sort_keys=True, separators=(',', ':'))))

        #main_table_id = self.metagenomic.insert_main_table('env_vif', mongo_data)

        main_table_id = self.meta.insert_none_table('sg_env_vif')

        update_info = {str(main_table_id): 'sg_env_vif'}

        #[geneset_table, gene_list] = metagenomic.export_geneset_table(data.geneset_id, data.method)

        options = {
            'otu_id': data.otu_id,
            'abund_file': data.otu_id,
            #'method': data.method,
            'update_info': json.dumps(update_info),
            'params': json.dumps(params_json, sort_keys=True, separators=(',', ':')),
            'group_id': data.group_id,
            #'geneset_id': data.geneset_id,
            'group_detail': data.group_detail,
            #'geneset_table': geneset_table,
            #'anno_type': data.anno_type,
            #'group': data.geneset_id,  # 为修改id转样品功能
            'env_labs': data.env_labs,
            'env_id': data.env_id,
            'env_file': data.env_id,
            'level':data.level_id,
            'viflim': int(data.viflim),
            'group_file': data.group_id
        }
        options['main_table_data'] = SON(mongo_data)

        to_file = []
        #to_file.append('metagenomic.export_group_table_by_detail(group)')
        to_file.append('meta.export_group_table_by_detail(group_file)')
        
        #to_file.append('metagenomic.export_float_env(env_file)')
        to_file.append('env.export_env_table(env_file)')
        #to_file.append('meta.export_otu_table_by_level(abund_file)')
        to_file.append('meta.export_otu_table_by_level(abund_file)')

        options['main_id'] = str(main_table_id)
        self.set_sheet_data(name=task_name, options=options,
                            main_table_name=main_table_name.strip().split("_")[0] + '/' + main_table_name,
                            module_type="workflow", params=params_json,
                            to_file=to_file)
        
        task_info = super(EnvVifAction, self).POST()
        if task_info['success']:
            task_info['content'] = {'ids': {'id': str(main_table_id), 'name': main_table_name}}
        return json.dumps(task_info)

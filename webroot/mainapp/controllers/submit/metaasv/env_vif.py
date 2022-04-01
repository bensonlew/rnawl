# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'

import web
import json
import datetime
from bson import ObjectId
from bson import SON
from mainapp.libs.param_pack import group_detail_sort
from mainapp.libs.signature import check_sig
from mainapp.controllers.project.metaasv_controller import MetaasvController


class EnvVifAction(MetaasvController):
    """
    Metaasv VIF方差膨胀因子分析
    """
    def __init__(self):
        super(EnvVifAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        default_argu = ['asv_id','level_id','submit_location', 'group_id', 'group_detail',"env_id", "env_labs", "viflim"]
        for argu in default_argu:
            if not hasattr(data, argu):
                info = {"success": False, "info": "parameters missing:%s" % argu}
                return json.dumps(info)

        task_name = 'metaasv.report.env_vif'
        module_type = 'workflow'
        task_type = 2

        otu_info = self.metaasv.get_otu_table_info(data.asv_id)
        task_info = self.metaasv.get_task_info(otu_info['task_id'])

        params_json = {
            'asv_id': data.asv_id,
            'level_id': int(data.level_id),
            'group_id': data.group_id,
            'group_detail': group_detail_sort(data.group_detail),
            'env_id': data.env_id,
            'env_labs': data.env_labs,
            'task_type': str(data.task_type),
            'submit_location': data.submit_location,
            'viflim': int(data.viflim)
        }
        group_id = data.group_id if data.group_id in ['all', 'All', 'ALL'] else ObjectId(data.group_id)
        mongo_data = [
            ('project_sn', task_info['project_sn']),
            ('task_id', task_info['task_id']),
            ('asv_id',ObjectId(data.asv_id)),
            ('status', 'start'),
            ('desc', '正在计算'),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        ]

        main_table_name = "Vif" + '_'  +  datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        mongo_data.append(('name', main_table_name))
        mongo_data.append(('params', json.dumps(params_json, sort_keys=True, separators=(',', ':'))))

        main_table_id = self.metaasv.insert_none_table('vif')
        update_info = {str(main_table_id): 'vif'}

        options = {
            'asv_id': data.asv_id,
            'abund_file': data.asv_id,
            'update_info': json.dumps(update_info),
            'params': json.dumps(params_json, sort_keys=True, separators=(',', ':')),
            'group_id': data.group_id,
            'group_detail': data.group_detail,
            'env_labs': data.env_labs,
            'env_file': data.env_id,
            'level':data.level_id,
            'viflim': int(data.viflim),
            'group_file': data.group_id
        }
        options['main_table_data'] = SON(mongo_data)

        to_file = []
        to_file.append('metaasv.export_group_table_by_detail(group_file)')
        to_file.append('metaasv_env.export_env_table(env_file)')
        to_file.append('metaasv.export_otu_table_by_level(abund_file)')

        options['main_id'] = str(main_table_id)
        self.set_sheet_data(name=task_name, options=options,
                            main_table_name=main_table_name.strip().split("_")[0] + '/' + main_table_name,
                            module_type=module_type, params=params_json,
                            to_file=to_file)
        
        task_info = super(EnvVifAction, self).POST()
        if task_info['success']:
            task_info['content'] = {'ids': {'id': str(main_table_id), 'name': main_table_name}}
        return json.dumps(task_info)

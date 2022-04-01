
# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'

import web
import json
import datetime
from bson import ObjectId
from mainapp.libs.param_pack import group_detail_sort
from mainapp.controllers.project.metaasv_controller import MetaasvController
from mainapp.libs.signature import check_sig


class VpaAction(MetaasvController):
    """
    metaasv vpa分析
    """
    def __init__(self):
        super(VpaAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        print(data)
        default_argu = [ "asv_id",'submit_location', 'group_detail','task_type',
                         'group_id', 'env_id','env_group_id', 'env_detail', 'level_id']
        for argu in default_argu:
            if not hasattr(data, argu):
                info = {"success": False, "info": "PARAMETERS MISSING: %s" , "variables":[ argu]}
                return json.dumps(info)

        task_name = 'metaasv.report.vpa'
        module_type = 'workflow'
        otu_info = self.metaasv.get_otu_table_info(data.asv_id)
        if not otu_info:
            info = {"success": False, "info": "OTU不存在，请确认参数是否正确！"}
            return json.dumps(info)
        task_info = self.metaasv.get_task_info(otu_info['task_id'])
        params_json = {
            'asv_id': data.asv_id,
            'level_id': int(data.level_id),
            'group_id': data.group_id,
            'group_detail': group_detail_sort(data.group_detail),
            'env_id': data.env_id,
            'env_detail': group_detail_sort(data.env_detail),
            'env_group_id': data.env_group_id,
            'submit_location': data.submit_location,
            'task_type': data.task_type
        }
        if hasattr(data, "env_labs"):
            params_json["env_labs"] = data.env_labs

        main_table_name = 'VPA_' + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        mongo_data = [
            ('project_sn', task_info['project_sn']),
            ('task_id', task_info['task_id']),
            ('asv_id', ObjectId(data.asv_id)),
            ('status', 'start'),
            ('desc', '正在计算'),
            ('name', main_table_name),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ("level_id", int(data.level_id)),
            ("params", json.dumps(params_json, sort_keys=True, separators=(',', ':')))
        ]
        main_table_id = self.metaasv.insert_main_table('vpa', mongo_data)
        update_info = {str(main_table_id): 'vpa'}
        options = {
            'otu_table': data.asv_id,
            'asv_id': data.asv_id,
            'level': int(data.level_id),
            'env_file': data.env_id,
            'env_labs': ','.join([','.join(i) for i in json.loads(data.env_detail).values()]),
            'env_detail' : data.env_detail,
            'env_group':data.env_detail,
            'group_detail': data.group_detail,
            'update_info': json.dumps(update_info),
            'main_id': str(main_table_id),
            'group_id': data.group_id,
            'group': data.group_id
        }

        to_file = ["metaasv.export_otu_table_by_level(otu_table)",
                    "metaasv.export_group_table_by_detail(group)", "metaasv_env.export_float_env_regression(env_file)"]
        to_file.append('metaasv.export_env_table_by_detail(env_group)')

        self.set_sheet_data(name=task_name, options=options, main_table_name="VPA/"+main_table_name,
                            module_type=module_type, to_file=to_file)
        task_info = super(VpaAction, self).POST()
        task_info['content'] = {
            'ids': {'id': str(main_table_id), 'name': main_table_name}}
        return json.dumps(task_info)
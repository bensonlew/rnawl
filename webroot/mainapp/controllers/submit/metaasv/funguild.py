# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'

import web
import json
from mainapp.libs.param_pack import group_detail_sort
from mainapp.controllers.project.metaasv_controller import MetaasvController
from mainapp.libs.signature import check_sig
from bson import ObjectId
import datetime


class FunguildAction(MetaasvController):
    """
    Metaasv funfuild分析
    """
    def __init__(self):
        super(FunguildAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        print(data)
        default_argu = ['asv_id', 'submit_location', 'group_detail', 'group_id']

        for argu in default_argu:
            if not hasattr(data, argu):
                info = {'success': False, 'info': 'parameters missing:%s' , "variables":[ argu]}
                return json.dumps(info)

        task_name = 'metaasv.report.funguild'
        task_type = 'workflow'
        otu_info = self.metaasv.get_otu_table_info(data.asv_id)
        if not otu_info:
            info = {"success": False, "info": "ASV不存在，请确认参数是否正确！"}
            return json.dumps(info)
        task_info = self.metaasv.get_task_info(otu_info['task_id'])
        params_json = {
            'asv_id': data.asv_id,
            'group_id': data.group_id,
            'group_detail': group_detail_sort(data.group_detail),
            'submit_location': data.submit_location,
            'task_type': data.task_type
        }
        if hasattr(data, "method"):
            params_json["method"] = data.method
        if hasattr(data, "others"):
            params_json["others"] = data.others

        main_table_name = 'FUNGuild_' + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]

        mongo_data = [
            ('project_sn', task_info['project_sn']),
            ('task_id', task_info['task_id']),
            ('asv_id', ObjectId(data.asv_id)),
            ('status', 'start'),
            ('desc', '正在计算'),
            ('name', main_table_name),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ("params", json.dumps(params_json, sort_keys=True, separators=(',', ':')))
        ]
        main_table_id = self.metaasv.insert_main_table('funguild', mongo_data)
        update_info = {str(main_table_id): 'funguild'}
        options = {
            'otu_table': data.asv_id,
            'level' : 9,
            'group_table': data.group_id,
            'group_detail': data.group_detail,
            'update_info': json.dumps(update_info),
            'main_id': str(main_table_id)
        }

        if hasattr(data, "method"):
            if data.method in ["none"]:
                options["method"] = ""
            else:
                options["method"] = data.method
        if hasattr(data, "others") :
            options["others"] = data.others

        to_file = ["metaasv.export_otu_table_by_level(otu_table)","metaasv.export_group_table_by_detail(group_table)"]

        self.set_sheet_data(name=task_name, options=options, main_table_name="FUNGuild/"+main_table_name,
                            module_type=task_type, to_file=to_file)
        task_info = super(FunguildAction, self).POST()
        task_info['content'] = {
            'ids': {'id': str(main_table_id), 'name': main_table_name}}
        return json.dumps(task_info)

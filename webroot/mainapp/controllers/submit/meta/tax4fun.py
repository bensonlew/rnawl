#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# @Author  : houshuang 2019/9/27

import web
import json
from bson.objectid import ObjectId
from mainapp.libs.param_pack import group_detail_sort, param_pack
from mainapp.models.mongo.meta import Meta
import datetime
from mainapp.controllers.project.meta_controller import MetaController
from mainapp.libs.signature import check_sig
from bson import SON


class Tax4funAction(MetaController):
    def __init__(self):
        super(Tax4funAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        print(data)
        default_argu = ['otu_id', 'submit_location', 'group_id', 'group_detail', 'method']
        for argu in default_argu:
            if not hasattr(data, argu):
                info = {'success': False, 'info': 'PARAMETERS MISSING: %s' , "variables":[ argu], "code" : "C2204601"}
                return json.dumps(info)
        if data.method not in ["", "sum", "average", "middle"]:
            info = {"success": False, "info": "对分组样本计算方式:%s错误!" , "variables":[ data.method], "code" : "C2204602"}
            return json.dumps(info)
        task_name = 'meta.report.tax4fun'
        task_type = 'workflow'
        otu_info = self.meta.get_otu_table_info(data.otu_id)
        if not otu_info:
            info = {"success": False, "info": "PARAMETERS ERROR: wrong value of otu_id, not found!", "code" : "C2204603"}
            return json.dumps(info)
        task_info = self.meta.get_task_info(otu_info['task_id'])
        params_json = {
            'otu_id': data.otu_id,
            'submit_location': data.submit_location,
            'group_id': data.group_id,
            'group_detail': group_detail_sort(data.group_detail),
            'method': data.method,
            'task_type': "reportTask"
        }
        main_table_name = 'Tax4fun_' + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        mongo_data = [
            ('project_sn', task_info['project_sn']),
            ('submit_location', data.submit_location),
            ('task_id', task_info['task_id']),
            ('otu_id', ObjectId(data.otu_id)),
            ('status', 'start'),
            ('desc', 'processing'),
            ('name', main_table_name),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        ]
        mongo_data.append(('params', json.dumps(params_json, sort_keys=True, separators=(',', ':'))))
        main_table_id = self.meta.insert_main_table('sg_tax4fun', mongo_data)
        update_info = {str(main_table_id): 'sg_tax4fun'}
        options = {
            'input_otu_id': data.otu_id,
            'in_otu_table': data.otu_id,
            'method': data.method,
            'group_detail': data.group_detail,
            'main_id': str(main_table_id),
            'task_id': task_info['task_id'],
            "level": int(9),
            'update_info': json.dumps(update_info)
        }
        to_file = 'meta.export_otu_table_by_detail(in_otu_table)' ###fix by qingchen.zhang@20191125
        self.set_sheet_data(name=task_name, options=options,
                            main_table_name=main_table_name.strip().split("_")[0] + '/' + main_table_name,
                            params=params_json,
                            module_type=task_type, to_file=to_file)
        task_info = super(Tax4funAction, self).POST()
        task_info['content'] = {'ids': {'id': str(main_table_id), 'name': main_table_name}}
        return json.dumps(task_info)

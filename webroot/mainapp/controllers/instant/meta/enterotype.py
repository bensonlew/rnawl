# -*- coding: utf-8 -*-
# __author__ = 'zouxuan'
import web
import json
from bson.objectid import ObjectId
from bson.errors import InvalidId
from mainapp.libs.param_pack import group_detail_sort
import types
from mainapp.models.mongo.meta import Meta
import datetime
from mainapp.controllers.project.meta_controller import MetaController
from bson import SON
from mainapp.libs.signature import check_sig
import time


class EnterotypeAction(MetaController):
    def __init__(self):
        super(EnterotypeAction, self).__init__(instant=True)

    @check_sig
    def POST(self):
        data = web.input()
        # default_argu = ['analysis_type', 'level_id', 'submit_location', 'group_id', 'group_detail','anno_id','geneset_id','second_level','anno_type']
        default_argu = ['submit_location', 'group_id', 'group_detail', 'otu_id', 'level_id',
                        'distance_method', 'task_type']
        for argu in default_argu:
            if not hasattr(data, argu):
                info = {'success': False, 'info': 'parameters missing:%s' % argu}
                return json.dumps(info)
        task_name = 'meta.report.enterotyping'
        module_type = 'workflow'
        main_table_name = 'Enterotypes_' + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        meta = Meta()
        otu_info = meta.get_otu_table_info(data.otu_id)
        if not otu_info:
            info = {"success": False, "info": "OTU不存在，请确认参数是否正确！", 'code':'C2200901'}
            return json.dumps(info)
        task_info = meta.get_task_info(otu_info['task_id'])
        group_detail = group_detail_sort(data.group_detail)
        params_json = {
            'otu_id': data.otu_id,
            'level_id': int(data.level_id),
            'submit_location': data.submit_location,
            'task_type': data.task_type,
            'group_id': data.group_id,
            'group_detail': group_detail,
            'distance_method': data.distance_method
        }
        # params_abu_json = {
        #     'geneset_id': data.geneset_id,
        #     'anno_type': data.anno_type,
        #     'method': data.method
        # }
        # abu_mongo = [
        #     ('project_sn', task_info['project_sn']),
        #     ('task_id', task_info['task_id']),
        #     ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        # ]
        # group_id = data.group_id if data.group_id in ['all', 'All', 'ALL'] else ObjectId(data.group_id)
        mongo_data = [
            ('project_sn', task_info['project_sn']),
            ('task_id', task_info['task_id']),
            ('status', 'start'),
            ('desc', '正在计算'),
            ('name', main_table_name),
            ('otu_id', ObjectId(data.otu_id)),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        ]
        mongo_data.append(('params', json.dumps(params_json, sort_keys=True, separators=(',', ':'))))
        main_table_id = meta.insert_none_table('sg_enterotype')
        update_info = {str(main_table_id): 'sg_enterotype'}
        options = {
            'otu_file': data.otu_id,
            'otu_id': data.otu_id,
            'distance_method': data.distance_method,
            'params': json.dumps(params_json, sort_keys=True, separators=(',', ':')),
            'group_id': data.group_id,
            'task_type': data.task_type,
            'group_detail': data.group_detail,
            'update_info':json.dumps(update_info),
            'group_table':data.group_id,
            'level':int(data.level_id),
            'main_id': str(main_table_id),
            'submit_location': data.submit_location,
            'main_table_data': SON(mongo_data)
        }
        to_file = ['meta.export_otu_table_by_detail(otu_file)',
                   'meta.export_group_table_by_detail(group_table)']
        self.set_sheet_data(name=task_name, options=options,
                            main_table_name=main_table_name.strip().split("_")[0] + '/' + main_table_name,
                            module_type=module_type, to_file=to_file)
        task_info = super(EnterotypeAction, self).POST()
        task_info['content'] = {'ids': {'id': str(main_table_id), 'name': main_table_name}}
        return json.dumps(task_info)

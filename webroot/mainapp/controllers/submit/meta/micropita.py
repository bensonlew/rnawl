# -*- coding: utf-8 -*-
# __author__ = ‘zhujuan'

import web
import json
from mainapp.controllers.project.meta_controller import MetaController
from mainapp.libs.signature import check_sig
from mainapp.libs.param_pack import group_detail_sort
from bson import ObjectId
from bson import SON
import datetime


class MicropitaAction(MetaController):
    def __init__(self):
        super(MicropitaAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        default_argu = ['otu_id', 'level_id', 'submit_location', 'group_detail', 'group_id', 'diversity_index',
                        'distance_method', 'specimen_n', 'task_type']
        for argu in default_argu:
            if not hasattr(data, argu):
                info = {'success': False, 'info': 'parameters missing:%s' % argu,  "code" : "C2202103"}
                return json.dumps(info)
        nu = 0
        group = []  # group = 1时表示，组内的样品个数都大于需要挑选的样品量，可以进行discriminant和 distinct的筛选
        group_detail = json.loads(data.group_detail)
        for i in group_detail.values():
            nu += len(i)
            if int(len(i)) > int(data.specimen_n):
                group.append(1)
            else:
                group.append(0)
        if int(nu) < int(data.specimen_n):
            info = {"success": False, "info": "丰度表的样本数目少于需要挑选的样品量，不可进行Micropita分析！", 'code':'C2202101'}
            return json.dumps(info)
        task_name = 'meta.report.micropita'
        module_type = 'workflow'
        otu_info = self.meta.get_otu_table_info(data.otu_id)
        if not otu_info:
            info = {"success": False, "info": "OTU不存在，请确认参数是否正确！", 'code':'C2202102'}
            return json.dumps(info)
        task_info = self.meta.get_task_info(otu_info['task_id'])
        params_json = {
            'otu_id': data.otu_id,
            'level_id': int(data.level_id),
            'group_id': data.group_id,
            'group_detail': group_detail_sort(data.group_detail),
            'diversity_index': data.diversity_index,
            'distance_method': data.distance_method,
            'specimen_n': data.specimen_n,
            'submit_location': data.submit_location,
            'task_type': data.task_type
        }
        if hasattr(data, 'target_feature'):
            params_json['target_feature'] = data.target_feature

        main_table_name = 'MicroPITA' + \
                          '_' + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        mongo_data = [
            ('project_sn', task_info['project_sn']),
            ('task_id', task_info['task_id']),
            ('otu_id', ObjectId(data.otu_id)),
            ('status', 'start'),
            ('desc', '正在计算'),
            ('name', main_table_name),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ("level_id", int(data.level_id)),
            ("params", json.dumps(params_json, sort_keys=True, separators=(',', ':')))
        ]
        main_table_id = self.meta.insert_none_table('sg_micropita')
        update_info = {str(main_table_id): 'sg_micropita'}
        options = {
            'otu_file': data.otu_id,
            'otu_id': data.otu_id,
            'level': int(data.level_id),
            'diversity_index': data.diversity_index,
            'distance_method': data.distance_method,
            'filter_nu': data.specimen_n,
            'group_id': data.group_id,
            'group_detail': data.group_detail,
            'update_info': json.dumps(update_info),
            "main_table_data": SON(mongo_data),
            "main_id": str(main_table_id)
        }
        if hasattr(data, 'target_feature'):
            options['target_feature'] = data.target_feature
        if 0 not in group and len(group_detail.keys()) > 1:
            options['group_file'] = data.group_id
            to_file = ["meta.export_otu_table_by_detail(otu_file)", "meta.export_group_table_by_detail(group_file)"]
        else:
            to_file = ["meta.export_otu_table_by_detail(otu_file)"]
        self.set_sheet_data(name=task_name, options=options, main_table_name="MicroPITA/" + main_table_name,
                            module_type=module_type, to_file=to_file)
        task_info = super(MicropitaAction, self).POST()
        if task_info['success']:
            task_info['content'] = {'ids': {'id': str(main_table_id), 'name': main_table_name}}
        return json.dumps(task_info)

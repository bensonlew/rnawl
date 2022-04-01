# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'

import web
import json
import datetime
from bson import ObjectId
from bson import SON
from mainapp.libs.param_pack import group_detail_sort,param_pack
from mainapp.libs.signature import check_sig
from mainapp.controllers.project.metaasv_controller import MetaasvController


class CompositionAction(MetaasvController):
    """
    metaasv 组成分析的circos图、Barpie图
    """
    def __init__(self):
        super(CompositionAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        default_argu = ['asv_id', 'submit_location', 'group_id', 'group_detail', 'group_method', 'graphic_type']
        for argu in default_argu:
            if not hasattr(data, argu):
                info = {"success": False, "info": "parameters missing:%s" % argu}
                return json.dumps(info)
        if data.group_method not in ["none", "sum", "average", "middle"]:
            variables = []
            variables.append(data.group_method)
            info = {"success": False, "info": "对分组样本计算方式:%s错误!" % data.group_method}
            return json.dumps(info)
        if not hasattr(data, 'level_id'):
            info = {"success": False, "info": "缺少level_id参数!"}
            return json.dumps(info)

        task_name = 'metaasv.report.composition'
        module_type = 'workflow'
        otu_info = self.metaasv.get_otu_table_info(data.asv_id)
        if not otu_info:
            info = {"success": False, "info": "OTU不存在，请确认参数是否正确！"}
            return json.dumps(info)
        task_info = self.metaasv.get_task_info(otu_info['task_id'])

        submit_location = data.submit_location
        if data.graphic_type in ['ternary','Ternary']:
            submit_location = 'ternary'
        params_json = {
            'group_id': data.group_id,
            'group_detail': group_detail_sort(data.group_detail),
            'group_method': data.group_method,
            'graphic_type': data.graphic_type,
            'task_type': str(data.task_type),
            'submit_location' : submit_location,
            'combine_value' : float(data.combine_value),
            'level_id' : int(data.level_id),
            'asv_id':data.asv_id
        }
        if hasattr(data, 'task_id'):
            params_json["task_id"] = data.task_id
        group_id = data.group_id if data.group_id in ['all', 'All', 'ALL'] else ObjectId(data.group_id)
        mongo_data = [
            ('asv_id',ObjectId(data.asv_id)),
            ('newick_id',None),
            ('project_sn', task_info['project_sn']),
            ('task_id', task_info['task_id']),
            ('status', 'start'),
            ('desc', '正在计算'),
            ('graphic_type', data.graphic_type),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
        ]
        if data.graphic_type in ["barpie"]:
            main_table_name = "CommunityBarPie" + '_' + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        elif data.graphic_type in ["circos"]:
            main_table_name = data.graphic_type.capitalize() + '_' + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        mongo_data.append(('name', main_table_name))
        mongo_data.append(('params', json.dumps(params_json, sort_keys=True, separators=(',', ':'))))

        if data.graphic_type in ["barpie"]:
            main_table_id = self.metaasv.insert_none_table('barpie')
            update_info = {str(main_table_id): 'barpie'}
        elif data.graphic_type in ["circos"]:
            main_table_id = self.metaasv.insert_none_table('circos')
            update_info = {str(main_table_id): 'circos'}

        options = {
            "asv_id": data.asv_id,
            "in_otu_table": data.asv_id,
            "group_detail": data.group_detail,
            "level": str(data.level_id),
            'update_info': json.dumps(update_info),
            'main_id': str(main_table_id),
            'main_table_data': SON(mongo_data),
            'combine_value' : data.combine_value,
            'graphic_type' : data.graphic_type,
            'group': data.group_id
            }
        if hasattr(data, "group_method"):
            if data.group_method in ["none"]:
                options["method"] = ""
            else:
                options["method"] = data.group_method
        to_file = ["metaasv.export_otu_table_by_level(in_otu_table)", "metaasv.export_group_table_by_detail(group)"]

        self.set_sheet_data(name=task_name, options=options,
                            main_table_name=data.graphic_type.capitalize() + '/' + main_table_name,
                            params=params_json,
                            module_type='workflow', to_file=to_file)
        task_info = super(CompositionAction, self).POST()
        if task_info['success']:
            task_info['content'] = {'ids': {'id': str(main_table_id), 'name': main_table_name}}
        return json.dumps(task_info)

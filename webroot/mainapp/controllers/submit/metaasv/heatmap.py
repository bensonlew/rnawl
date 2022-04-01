# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'

import web
import json
import datetime
from mainapp.controllers.project.metaasv_controller import MetaasvController
from mainapp.libs.param_pack import group_detail_sort
from mainapp.libs.signature import check_sig
from bson import ObjectId
from bson import SON


class HeatmapAction(MetaasvController):
    """
    metaasv 群落heatmap分析
    """
    def __init__(self):
        super(HeatmapAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        postArgs = ["asv_id", "level_id", "group_id", "group_detail", "top", "species_method", "task_type", "sample_method", "group_method"]
        for arg in postArgs:
            if not hasattr(data, arg):
                info = {'success': False, 'info': '%sparameters missing!', "variables":[arg]}
                return json.dumps(info)
        task_name = 'metaasv.report.composition_heatmap'
        otu_info = self.metaasv.get_otu_table_info(data.asv_id)
        if not otu_info:
            info = {"success": False, "info": "OTU不存在，请确认参数是否正确！"}
            return json.dumps(info)
        if data.species_method not in ["average", "single", "complete", "", "none"]:
            variables = []
            variables.append(data.species_method)
            info = {'success': False, "info": "参数species_method的值为{}，应该为average，single或者complete".format(data.species_method)}
            return json.dumps(info)
        task_info = self.metaasv.get_task_info(otu_info['task_id'])
        params_json = {
            'submit_location': data.submit_location,
            'asv_id': data.asv_id,
            'level_id': int(data.level_id),
            'group_detail': group_detail_sort(data.group_detail),
            'group_id': data.group_id,
            'task_type': data.task_type,
            'top': data.top,
            "species_method": data.species_method,
            "sample_method": data.sample_method,
        }
        if hasattr(data, "color_level"):
            params_json["color_level"] = int(data.color_level)
        if hasattr(data, "group_method"):
            params_json["group_method"] = str(data.group_method)
        if hasattr(data, "add_Algorithm"):
            params_json["add_Algorithm"] = str(data.add_Algorithm)
        if hasattr(data, "sample_add_Algorithm"):
            params_json["sample_add_Algorithm"] = str(data.sample_add_Algorithm)
        main_table_name = 'CommunityHeatmap_' + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        mongo_data = [
            ('project_sn', task_info['project_sn']),
            ('task_id', task_info['task_id']),
            ('asv_id', ObjectId(data.asv_id)),
            ('name', main_table_name),
            ("params", json.dumps(params_json, sort_keys=True, separators=(',', ':'))),
            ('status', 'start'),
            ('desc', '正在计算'),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
        ]
        main_table_id = self.metaasv.insert_none_table('heatmap')
        update_info = {str(main_table_id): 'heatmap'}
        options = {
            "input_otu_id": data.asv_id,
            "in_otu_table": data.asv_id,
            "group_detail": data.group_detail,
            "level": str(data.level_id),
            'update_info': json.dumps(update_info),
            'main_id': str(main_table_id),
            "species_number": data.top,  # 筛选物种参数
            'main_table_data': SON(mongo_data),
            'group': data.group_id
        }
        if hasattr(data, "species_method"):
            if data.species_method in ["none"]:
                options["method"] = ""
            else:
                options["method"] = data.species_method
        if hasattr(data, "sample_method"):# 样本聚类方式
            if data.sample_method in ["none"]:
                options["sample_method"] = ""
            else:
                options["sample_method"] = data.sample_method
        if hasattr(data, "group_method"):
            if data.group_method in ["none"]:
                options["add_Algorithm"] = ""
            else:
                options["add_Algorithm"] = data.group_method
        if hasattr(data, "color_level"):
            options["level_color"] = str(data.color_level)
        if hasattr(data, "add_Algorithm"):
            options["species_distance"] = str(data.add_Algorithm)
        if hasattr(data, "sample_add_Algorithm"):
            options["sample_distance"] = str(data.sample_add_Algorithm)

        to_file = ["metaasv.export_otu_table_by_level(in_otu_table)", "metaasv.export_group_table_by_detail(group)"]
        self.set_sheet_data(name=task_name, options=options, main_table_name="CommunityAnalysis/" + main_table_name,
                            module_type='workflow', to_file=to_file)
        task_info = super(HeatmapAction, self).POST()
        if task_info['success']:
            task_info['content'] = {
                'ids': {
                    'id': str(main_table_id),
                    'name': main_table_name
                }}
        return json.dumps(task_info)

# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'
import web
import json
from mainapp.libs.param_pack import group_detail_sort
from mainapp.controllers.project.metaasv_controller import MetaasvController
from mainapp.libs.signature import check_sig
from bson import ObjectId
import datetime


class RegressionAction(MetaasvController):
    """
    metaasv 排序回归分析
    """
    def __init__(self):
        super(RegressionAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        default_argu = ['asv_id', 'level_id', 'submit_location',
                        'group_detail', 'group_id', 'env_id', 'env_labs']

        for argu in default_argu:
            if not hasattr(data, argu):
                info = {'success': False, 'info': 'parameters missing:%s' , "variables":[ argu]}
                return json.dumps(info)
        task_name = 'metaasv.report.regression'
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
            'env_labs': data.env_labs,
            'submit_location': data.submit_location,
            'task_type': str(data.task_type)
        }

        if hasattr(data, "diversity_type"):
            params_json['diversity_type'] = data.diversity_type

        if hasattr(data, "distance_type"):
            params_json['distance_type'] = data.distance_type

        if hasattr(data,"diversity_analysis_type"):
            params_json['diversity_analysis_type'] = data.diversity_analysis_type

        main_table_name = 'EnvironmentalRegression_' + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
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
        main_table_id = self.metaasv.insert_main_table('regression', mongo_data)
        update_info = {str(main_table_id): 'regression'}
        options = {
            'otu_table': data.asv_id,
            'asv_id': data.asv_id,
            'level': int(data.level_id),
            'envtable': data.env_id,
            'env_labs': data.env_labs,
            'group_table': data.group_id,
            'group_detail': data.group_detail,
            'update_info': json.dumps(update_info),
            'main_id': str(main_table_id)
        }
        if hasattr(data, "diversity_type"):
            options['diversity_type'] = data.diversity_type

        if hasattr(data, "distance_type"):
            options['distance_type'] = data.distance_type

        if hasattr(data,"diversity_analysis_type"):
            options['diversity_analysis_type'] = data.diversity_analysis_type

        to_file = ["metaasv.export_otu_table_by_level(otu_table)","metaasv.export_group_table_by_detail(group_table)",
                   "metaasv_env.export_float_env_regression(envtable)"]

        self.set_sheet_data(name=task_name, options=options, main_table_name="Environmental_Regression/"+main_table_name,
                            module_type=module_type, to_file=to_file)
        task_info = super(RegressionAction, self).POST()
        task_info['content'] = {
            'ids': {'id': str(main_table_id), 'name': main_table_name}}
        return json.dumps(task_info)

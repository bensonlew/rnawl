# -*- coding: utf-8 -*-
# __author__ = 'shenghe'
import web
import json
import datetime
from bson import SON
from mainapp.libs.param_pack import group_detail_sort
from mainapp.controllers.project.metaasv_controller import MetaasvController
from mainapp.models.mongo.metaasv import Metaasv
from mainapp.libs.signature import check_sig
from bson import ObjectId


class HclusterAction(MetaasvController):
    """
    Metaasv 样本层级分析
    """
    def __init__(self):
        super(HclusterAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        default_argu = ['asv_id', 'level_id', 'group_id', 'group_detail',
                        'distance_algorithm', 'hcluster_method', 'submit_location']
        for argu in default_argu:
            if not hasattr(data, argu):
                info = {'success': False, 'info': 'parameters missing!: %s' , "variables":[ argu]}
                return json.dumps(info)
        if hasattr(data, "others_value"):
            if float(data.others_value) < 0 or float(data.others_value) > 1:
                variables = [data.others_value]
                info = {'success': False, "info": "柱形图Others合并值为{}，大小范围应当为0-1".format(data.others_value), 'variables': variables,}
                return json.dumps(info)
        task_name = 'metaasv.report.hcluster'
        task_type = 'workflow'
        meta = Metaasv()
        otu_info = meta.get_otu_table_info(data.asv_id)
        if not otu_info:
            info = {"success": False, "info": "ASV不存在，请确认参数是否正确！"}
            return json.dumps(info)
        task_info = meta.get_task_info(otu_info['task_id'])
        params_json = {
            'asv_id': data.asv_id,
            'level_id': int(data.level_id),
            'group_id': data.group_id,
            'group_detail': group_detail_sort(data.group_detail),
            'distance_algorithm': data.distance_algorithm,
            'hcluster_method': data.hcluster_method,
            'submit_location': data.submit_location,
            'task_type': data.task_type
        }
        if hasattr(data, "others_value"):
            params_json["others_value"] = data.others_value
        main_table_name = 'Hcluster' + data.hcluster_method.capitalize() + \
            '_' + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
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
        main_table_id = meta.insert_none_table('hcluster')
        update_info = {str(main_table_id): 'hcluster'}
        options = {
            'otu_table': data.asv_id,
            'asv_id': data.asv_id,
            'level': int(data.level_id),
            'dist_method': data.distance_algorithm,
            'hcluster_method': data.hcluster_method,
            'group_detail': data.group_detail,
            'update_info': json.dumps(update_info),
            'params': json.dumps(params_json, sort_keys=True, separators=(',', ':')),
            'main_id': str(main_table_id),
            "group_file": data.group_id,
            "main_table_data": SON(mongo_data)
        }
        if hasattr(data, "others_value"):
            options["others_value"] = float(data.others_value)
        to_file = ['metaasv.export_otu_table_by_detail(otu_table)', 'metaasv.export_group_table_by_detail(group_file)']
        self.set_sheet_data(name=task_name, options=options, main_table_name="HclusterAnalysis/" + main_table_name,
                            module_type=task_type, to_file=to_file)
        task_info = super(HclusterAction, self).POST()
        task_info['content'] = {'ids': {'id': str(main_table_id), 'name': main_table_name}}
        return json.dumps(task_info)

# -*- coding: utf-8 -*-
# __author__ = 'shenghe'
import web
import json
import datetime
from bson import SON
from mainapp.libs.param_pack import group_detail_sort
from mainapp.controllers.project.meta_controller import MetaController
from mainapp.models.mongo.meta import Meta
from mainapp.libs.signature import check_sig
from bson import ObjectId


class HclusterAction(MetaController):

    def __init__(self):
        super(HclusterAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        default_argu = ['otu_id', 'level_id', 'group_id', 'group_detail',
                        'distance_algorithm', 'hcluster_method', 'submit_location']
        for argu in default_argu:
            if not hasattr(data, argu):
                info = {'success': False, 'info': 'parameters missing!: %s' , "variables":[ argu], "code" : "C2200202"}
                return json.dumps(info)
        # by houshuang 20191017 >>>
        if hasattr(data, "others_value"):
            if float(data.others_value) < 0 or float(data.others_value) > 1:
                variables = [data.others_value]
                info = {'success': False, "info": "柱形图Others合并值为{}，大小范围应当为0-1".format(data.others_value), 'variables': variables, "code" : "C2200203"}
                return json.dumps(info)
        # <<<
        task_name = 'meta.report.hcluster'
        task_type = 'workflow'
        meta = Meta()
        otu_info = meta.get_otu_table_info(data.otu_id)
        if not otu_info:
            info = {"success": False, "info": "OTU不存在，请确认参数是否正确！", 'code':'C2200201'}
            return json.dumps(info)
        task_info = meta.get_task_info(otu_info['task_id'])
        params_json = {
            'otu_id': data.otu_id,
            'level_id': int(data.level_id),
            'group_id': data.group_id,
            'group_detail': group_detail_sort(data.group_detail),
            'distance_algorithm': data.distance_algorithm,
            'hcluster_method': data.hcluster_method,
            'submit_location': data.submit_location,
            'task_type': data.task_type
        }
        # by houshuang 20190918 >>>
        if hasattr(data, "others_value"):
            params_json["others_value"] = data.others_value
        # <<<
        main_table_name = 'Hcluster' + data.hcluster_method.capitalize() + \
            '_' + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        mongo_data = [
            ('project_sn', task_info['project_sn']),
            ('task_id', task_info['task_id']),
            ('table_id', ObjectId(data.otu_id)),
            ('table_type', 'dist'),
            ('tree_type', 'cluster'),
            ('status', 'start'),
            ('desc', '正在计算'),
            ('name', main_table_name),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ("level_id", int(data.level_id)),
            ("params", json.dumps(params_json, sort_keys=True, separators=(',', ':')))
        ]
        main_table_id = meta.insert_none_table('sg_newick_tree')
        update_info = {str(main_table_id): 'sg_newick_tree'}
        options = {
            'otu_table': data.otu_id,
            'otu_id': data.otu_id,
            'level': int(data.level_id),
            'dist_method': data.distance_algorithm,
            'hcluster_method': data.hcluster_method,
            'group_detail': data.group_detail,
            'update_info': json.dumps(update_info),
            'params': json.dumps(params_json, sort_keys=True, separators=(',', ':')),
            'main_id': str(main_table_id),
            "main_table_data": SON(mongo_data)
        }
        # by houshuang 20190918 >>>
        if hasattr(data, "others_value"):
            options["others_value"] = float(data.others_value)
        # <<<
        to_file = 'meta.export_otu_table_by_detail(otu_table)'
        self.set_sheet_data(name=task_name, options=options, main_table_name="HclusterAnalysis/" + main_table_name,
                            module_type=task_type, to_file=to_file)
        task_info = super(HclusterAction, self).POST()
        task_info['content'] = {'ids': {'id': str(main_table_id), 'name': main_table_name}}
        return json.dumps(task_info)

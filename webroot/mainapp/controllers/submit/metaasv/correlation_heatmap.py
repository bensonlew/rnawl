# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'
import web
import json
from mainapp.controllers.project.metaasv_controller import MetaasvController
from mainapp.libs.param_pack import group_detail_sort
from mainapp.libs.signature import check_sig
from bson import ObjectId
from bson import SON
import datetime


class CorrelationHeatmapAction(MetaasvController):
    """
    metaasv 相关性heatmap图分析接口
    """
    def __init__(self):
        super(CorrelationHeatmapAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        default_argu = ['asv_id', 'level_id', 'submit_location', "group_id", "env_id", "env_labs", "top_species"]
        if not hasattr(data, 'env_id'):
            info = {'success': False, 'info': '缺少环境因子参数!'}
            return json.dumps(info)
        for argu in default_argu:
            if not hasattr(data, argu):
                info = {'success': False, 'info': '%s参数缺少!' % argu}
                return json.dumps(info)
        if len(str(data.env_labs).split(",")) < 2:
            info = {'success': False, 'info': '相关性Heatmap分析环境因子数量必须大于等于2!'}
            return json.dumps(info)
        task_name = 'metaasv.report.correlation_heatmap'
        module_type = 'workflow'
        otu_info = self.metaasv.get_otu_table_info(data.asv_id)
        if not otu_info:
            info = {"success": False, "info": "ASV不存在，请确认参数是否正确！!"}
            return json.dumps(info)
        task_info = self.metaasv.get_task_info(otu_info['task_id'])

        params_json = {
            "asv_id": data.asv_id,
            "level_id": int(data.level_id),
            "submit_location": data.submit_location,
            "task_type": str(data.task_type),
            "group_detail": group_detail_sort(data.group_detail),
            "group_id": data.group_id,
            "env_id": data.env_id,
            "env_labs": data.env_labs
        }
        if hasattr(data, "method"):
            params_json["method"] = data.method
        if hasattr(data, "sample_method"):
            params_json["sample_method"] = data.sample_method
        if hasattr(data, "species_method"):
            params_json["species_method"] = data.species_method
        if hasattr(data, "top_species"):
            params_json["top_species"] = int(data.top_species)
        if hasattr(data, "env_distance"):
            params_json["env_distance"] = data.env_distance
        if hasattr(data, "species_distance"):
            params_json["species_distance"] = data.species_distance

        level_name = ["Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "ASV"]
        main_table_name = 'CorrelationHeatmap' + level_name[int(data.level_id) - 1] + "_" + \
            datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        mongo_data = [
            ('project_sn', task_info['project_sn']),
            ('task_id', task_info['task_id']),
            ('asv_id', ObjectId(data.asv_id)),
            ("env_id", ObjectId(data.env_id)),
            ('status', 'start'),
            ('name', main_table_name),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ("level_id", int(data.level_id)),
            ("params", json.dumps(params_json, sort_keys=True, separators=(',', ':')))
        ]
        main_table_id = self.metaasv.insert_main_table('cor_heatmap', mongo_data)
        update_info = {str(main_table_id): 'cor_heatmap'}

        options = {
            "otu_file": data.asv_id,
            "level": int(data.level_id),
            "env_file": data.env_id,
            'update_info': json.dumps(update_info),
            "group_detail": data.group_detail,
            "main_id": str(main_table_id),
            "main_table_data": SON(mongo_data),
            # "env_distance": data.env_distance,
            # "species_distance": data.species_distance,
            # "env_cluster": data.env_cluster,
            # "species_cluster": data.species_cluster
        }
        del params_json["level_id"]
        del params_json["group_detail"]
        del params_json["submit_location"]
        del params_json["task_type"]
        if hasattr(data, "method"):
            if data.method == "pearson":
                params_json["method"] = "pearsonr"
            else:
                params_json["method"] = "spearmanr"
        if hasattr(data, "sample_method"):
            params_json["sample_method"] = data.species_method
            if params_json["sample_method"] == "none":
                params_json["env_cluster"] = ""
            else:
                params_json["env_cluster"] = data.sample_method
            del params_json["sample_method"]
        if hasattr(data, "species_method"):
            params_json["species_method"] = data.species_method
            if params_json["species_method"] == "none":
                params_json["species_cluster"] = ""
            else:
                params_json["species_cluster"] = data.species_method
            del params_json["species_method"]
        options.update(params_json)
        print(options)
        to_file = ['metaasv.export_otu_table_by_detail(otu_file)', "metaasv_env.export_float_env(env_file)"]
        self.set_sheet_data(name=task_name, options=options, main_table_name="CorrelationHeatmap/" + main_table_name,
                            module_type=module_type, to_file=to_file)
        task_info = super(CorrelationHeatmapAction, self).POST()
        if task_info['success']:
            task_info['content'] = {
                'ids': {
                    'id': str(main_table_id),
                    'name': main_table_name
                }}
        return json.dumps(task_info)

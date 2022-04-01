# -*- coding: utf-8 -*-
# __author__ = 'qindanhua'
import web
import json
from mainapp.controllers.project.meta_controller import MetaController
from mainapp.libs.param_pack import group_detail_sort
from mainapp.libs.signature import check_sig
from bson import ObjectId
from bson import SON
import datetime


class PearsonCorrelationAction(MetaController):
    """
    pearson 相关系数分析接口
    """

    def __init__(self):
        super(PearsonCorrelationAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        default_argu = ['otu_id', 'level_id', 'submit_location', "group_id", "env_id", "env_labs", "top_species"]
        if not hasattr(data, 'env_id'):  # modified by hongdongxuan 20170310
            info = {'success': False, 'info': '缺少环境因子参数!'}
            return json.dumps(info)
        for argu in default_argu:
            if not hasattr(data, argu):
                info = {'success': False, 'info': '%s参数缺少!' % argu}
                return json.dumps(info)
        if len(str(data.env_labs).split(",")) < 2:
            info = {'success': False, 'info': '相关性Heatmap分析环境因子数量必须大于等于2!'}
            return json.dumps(info)
        task_name = 'meta.report.pearson_correlation'
        task_type = 'workflow'
        otu_info = self.meta.get_otu_table_info(data.otu_id)
        if not otu_info:
            info = {"success": False, "info": "OTU不存在，请确认参数是否正确！!"}
            return json.dumps(info)
        task_info = self.meta.get_task_info(otu_info['task_id'])

        # print(data.top_species)
        params_json = {
            "otu_id": data.otu_id,
            "level_id": int(data.level_id),
            "submit_location": data.submit_location,
            "task_type": data.task_type,
            "group_detail": group_detail_sort(data.group_detail),
            "group_id": data.group_id,
            "env_id": data.env_id,
            "env_labs": data.env_labs,
            # "method": "pearsonr"
        }
        method_name = "Pearson"
        if hasattr(data, "method"):
            # print(data.method)
            params_json["method"] = data.method
            if data.method == "pearsonr":
                method_name = "Pearson"
            else:
                method_name = "Spearman"
        if hasattr(data, "env_cluster"):
            params_json["env_cluster"] = data.env_cluster
        if hasattr(data, "species_cluster"):
            params_json["species_cluster"] = data.species_cluster
        if hasattr(data, "top_species"):
            params_json["top_species"] = data.top_species
        # self.options["params"] = str(self.options)
        level_name = ["Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "OTU"]  # add by hongdongxuan 20170322
        main_table_name = method_name + 'Correlation' + level_name[int(data.level_id) - 1] + "_" + \
            datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        mongo_data = [
            ('project_sn', task_info['project_sn']),
            ('task_id', task_info['task_id']),
            ('otu_id', ObjectId(data.otu_id)),
            ("env_id", ObjectId(data.env_id)),
            ('status', 'start'),
            ('name', main_table_name),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ("level_id", int(data.level_id)),
            ("params", json.dumps(params_json, sort_keys=True, separators=(',', ':')))
        ]
        main_table_id = self.meta.insert_none_table('sg_species_env_correlation')
        update_info = {str(main_table_id): 'sg_species_env_correlation'}

        options = {
            "otu_file": data.otu_id,
            "level": int(data.level_id),
            "env_file": data.env_id,
            'update_info': json.dumps(update_info),
            "group_detail": data.group_detail,
            "corr_id": str(main_table_id),
            "main_table_data": SON(mongo_data)
        }
        del params_json["level_id"]
        del params_json["group_detail"]
        options.update(params_json)
        # print("lllllloptionlllll")
        # print(options)
        # print("ooooooooptionoooooo")
        to_file = ['meta.export_otu_table_by_detail(otu_file)', "env.export_float_env(env_file)"]
        self.set_sheet_data(name=task_name, options=options, main_table_name="CorrelationHeatmap/" + main_table_name,
                            module_type=task_type, to_file=to_file)
        task_info = super(PearsonCorrelationAction, self).POST()
        if task_info['success']:
            task_info['content'] = {
                'ids': {
                    'id': str(main_table_id),
                    'name': main_table_name
                }}
        return json.dumps(task_info)

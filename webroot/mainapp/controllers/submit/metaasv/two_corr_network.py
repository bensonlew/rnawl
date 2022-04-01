# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'

import web
import json
from bson.objectid import ObjectId
from mainapp.libs.param_pack import group_detail_sort
import datetime
from mainapp.controllers.project.metaasv_controller import MetaasvController
from mainapp.libs.signature import check_sig


class TwoCorrNetworkAction(MetaasvController):
    """
    Metaasv 双因素相关性分析
    """
    def __init__(self):
        super(TwoCorrNetworkAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        print data
        default_argu = ['submit_location', 'group_id', 'group_detail',  'ratio_method','color_level', 'task_type', 'level_id', 'asv_id']
        for argu in default_argu:
            if not hasattr(data, argu):
                info = {"success": False, "info": "PARAMETERS MISSING: %s" , "variables":[ argu],}
                return json.dumps(info)

        task_name = "metaasv.report.network_corfd"
        module_type = "workflow"
        task_type = data.task_type

        otu_info = self.metaasv.get_otu_table_info(data.asv_id)
        if not otu_info:
            info = {"success": False, "info": "ASV不存在，请确认参数是否正确！"}
            return json.dumps(info)
        task_info = self.metaasv.get_task_info(otu_info['task_id'])


        # if float(data.cor_value) > 1 or float(data.cor_value) < 0:
        #     info = {"success": False, "info": "参数错误：相关系数阈值范围应该在[0-1]", 'variables': ''}
        #     return json.dumps(info)
        # if float(data.significance) > 1 or float(data.significance) <= 0:
        #     info = {"success": False, "info": "参数错误：P-value范围应该在（0,1]", 'variables': ''}
        #     return json.dumps(info)

        group_detail = group_detail_sort(data.group_detail)
        params_json = {
            "submit_location": data.submit_location,
            "task_type": task_type,
            "group_id": data.group_id,
            "group_detail": group_detail,
            "ratio_method": data.ratio_method,
            "level_id" : int(data.level_id),
            'asv_id': data.asv_id,
            'env_id': data.env_id,
            'env_labs': data.env_labs,
            'color_level' : data.color_level
        }
        if hasattr(data, "cor_value"):
            params_json["cor_value"] = data.cor_value
        if hasattr(data, "significance"):
            params_json["significance"] = float(data.significance)
        if hasattr(data, "significance"):
            params_json["significance"] = float(data.significance)
        if hasattr(data, 'abundance'):
            params_json['abundance'] = int(data.abundance)
        name = "TwowayCorrNetwork_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        task_id = task_info['task_id']
        project_sn = task_info['project_sn']

        mongo_data = [
            ('project_sn', project_sn),
            ('task_id',task_id),
            ("asv_id", ObjectId(data.asv_id)),
            ("status", "start"),
            ("name", name),
            ("desc", "processing"),
            ("created_ts", datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
        ]

        options = {
            "otu_table" : data.asv_id ,
            "level":  int(data.level_id),
            "asv_id" : data.asv_id,
            "env_id": data.env_id,
            "env_file" : data.env_id,
            "group_table" : data.group_id,
            "group_detail" :  data.group_detail,
            "env_labs" : data.env_labs,
            "coefficient": data.ratio_method,
            "color_level" : data.color_level
        }
        if hasattr(data, "significance"):
            options["pvalue"] = float(data.significance)
        if hasattr(data, "cor_value"):
            options["coefficient_value"] = float(data.cor_value)
        if hasattr(data, "abundance"):
            options["fac1_top"] = int(data.abundance)

        to_file = ["metaasv.export_otu_table_by_level(otu_table)" ,
                   "metaasv_env.export_float_env_regression(env_file)", "metaasv.export_group_table_by_detail(group_table)" ]


        mongo_data.append(('params', json.dumps(params_json, sort_keys=True, separators=(",", ":"))))
        main_table_id = self.metaasv.insert_main_table("two_corr_network", mongo_data)

        update_info = {str(main_table_id): "two_corr_network"}
        options["update_info"] = json.dumps(update_info)
        options["main_table_id"] = str(main_table_id)
        self.set_sheet_data(name=task_name, options=options, main_table_name="TwowayCorrNetwork/" + name,
                            module_type=module_type, to_file=to_file )  #task_id=task_id, params=params_json

        task_info = super(TwoCorrNetworkAction, self).POST()
        if task_info['success']:
            task_info['content'] = {'ids': {'id': str(main_table_id), 'name': name}}
        return json.dumps(task_info)


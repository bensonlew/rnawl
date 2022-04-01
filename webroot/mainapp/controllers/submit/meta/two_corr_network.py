# -*- coding: utf-8 -*-


import web
import json
from bson.objectid import ObjectId

from mainapp.libs.param_pack import group_detail_sort
import types
import datetime
from mainapp.controllers.project.meta_controller import MetaController
from bson import SON
import os
from mainapp.libs.signature import check_sig



class TwoCorrNetworkAction(MetaController):
    def __init__(self):
        super(TwoCorrNetworkAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        print data
        default_argu = ['submit_location', 'group_id', 'group_detail',  'ratio_method','color_level',
                        'cor_value', 'significance', 'task_type', 'level_id', 'otu_id', 'abundance']  #'geneset_id', 'method', ,"factor1", "factor2"
        for argu in default_argu:
            if not hasattr(data, argu):
                info = {"success": False, "info": "PARAMETERS MISSING: %s" , "variables":[ argu], "code" : "C2204301"}
                return json.dumps(info)

        task_name = "meta.report.network_corfd"  # 调用workflow
        module_type = "workflow"
        task_type = data.task_type

        otu_info = self.meta.get_otu_table_info(data.otu_id)
        if not otu_info:
            info = {"success": False, "info": "OTU不存在，请确认参数是否正确！", "code" : "C2204302"}
            return json.dumps(info)
        task_info = self.meta.get_task_info(otu_info['task_id'])


        if float(data.cor_value) > 1 or float(data.cor_value) < 0:
            info = {"success": False, "info": "参数错误：相关系数阈值范围应该在[0-1]", 'variables': '', "code" : "C2204303"}
            return json.dumps(info)
        if float(data.significance) > 1 or float(data.significance) <= 0:
            info = {"success": False, "info": "参数错误：P-value范围应该在（0,1]", 'variables': '', "code" : "C2204304"}
            return json.dumps(info)


        group_detail = group_detail_sort(data.group_detail)
        #group_id = "all" if data.group_id in ['all', 'All', 'ALL'] else ObjectId(data.group_id)
        params_json = {
            "submit_location": data.submit_location,
            "task_type": task_type,
            "group_id": data.group_id,
            "group_detail": group_detail,
            "ratio_method": data.ratio_method,
            "cor_value": data.cor_value,
            "significance": float(data.significance),
            "level_id" : data.level_id,
            'otu_id': data.otu_id,
            'env_id': data.env_id,
            'env_labs': data.env_labs,
            'abundance' : data.abundance,
            'color_level' : data.color_level
        }

        name = "TwowayCorrNetwork_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        task_id = task_info['task_id']
        project_sn = task_info['project_sn']

        mongo_data = [
            ('project_sn', project_sn),
            ('task_id',task_id),
            ("otu_id", ObjectId(data.otu_id)),
            ("status", "start"),
            ("name", name),
            ("desc", "processing"),
            ("created_ts", datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            #("specimen", ",".join(samples_ids)),
        ]

        options = {
            "otu_table" : data.otu_id ,
            "level":  int(data.level_id),
            "otu_id" : data.otu_id,
            "env_id": data.env_id,
            "env_file" : data.env_id,
            "group_table" : data.group_id,
            "group_detail" :  data.group_detail,
            "env_labs" : data.env_labs,
            "fac1_top" : data.abundance,
            "coefficient": data.ratio_method,
            "coefficient_value": data.cor_value,
            "pvalue": data.significance,
            "color_level" : data.color_level
        }

        to_file = ["meta.export_otu_table_by_level(otu_table)" ,
                   "env.export_float_env_regression(env_file)", "meta.export_group_table_by_detail(group_table)" ]


        mongo_data.append(('params', json.dumps(params_json, sort_keys=True, separators=(",", ":"))))
        main_table_id = self.meta.insert_main_table("sg_two_corr_network", mongo_data)

        update_info = {str(main_table_id): "sg_two_corr_network"}
        options["update_info"] = json.dumps(update_info)
        options["main_table_id"] = str(main_table_id)
        self.set_sheet_data(name=task_name, options=options, main_table_name="TwowayCorrNetwork/" + name,
                            module_type=module_type, to_file=to_file )  #task_id=task_id, params=params_json

        task_info = super(TwoCorrNetworkAction, self).POST()
        if task_info['success']:
            task_info['content'] = {'ids': {'id': str(main_table_id), 'name': name}}
        return json.dumps(task_info)


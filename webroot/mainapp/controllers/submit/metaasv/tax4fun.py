#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# __author__ = "qingchen.zhang"

import web
import json
from bson.objectid import ObjectId
from mainapp.libs.param_pack import group_detail_sort, param_pack
import datetime
from mainapp.controllers.project.metaasv_controller import MetaasvController
from mainapp.libs.signature import check_sig


class Tax4funAction(MetaasvController):
    def __init__(self):
        """
        Metaasv Tax4Fun分析
        """
        super(Tax4funAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        print(data)
        default_argu = ["asv_id", "submit_location", "group_id", "group_detail", "method"]
        for argu in default_argu:
            if not hasattr(data, argu):
                info = {"success": False, "info": "PARAMETERS MISSING: %s" , "variables":[ argu]}
                return json.dumps(info)
        if data.method not in ["", "sum", "average", "middle"]:
            info = {"success": False, "info": "对分组样本计算方式:%s错误!" , "variables":[ data.method]}
            return json.dumps(info)
        task_name = "metaasv.report.tax4fun"
        module_type = "workflow"
        otu_info = self.metaasv.get_otu_table_info(data.asv_id)
        if not otu_info:
            info = {"success": False, "info": "PARAMETERS ERROR: wrong value of otu_id, not found!", "code" : "C2204603"}
            return json.dumps(info)
        task_info = self.metaasv.get_task_info(otu_info["task_id"])
        params_json = {
            "asv_id": data.asv_id,
            "submit_location": data.submit_location,
            "group_id": data.group_id,
            "group_detail": group_detail_sort(data.group_detail),
            "method": data.method,
            "task_type": str("2")
        }
        main_table_name = "Tax4fun_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        mongo_data = [
            ("project_sn", task_info["project_sn"]),
            ("submit_location", data.submit_location),
            ("task_id", task_info["task_id"]),
            ("asv_id", ObjectId(data.asv_id)),
            ("status", "start"),
            ("desc", "processing"),
            ("name", main_table_name),
            ("created_ts", datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        ]
        mongo_data.append(("params", json.dumps(params_json, sort_keys=True, separators=(",", ":"))))
        main_table_id = self.metaasv.insert_main_table("tax4fun", mongo_data)
        update_info = {str(main_table_id): "tax4fun"}
        options = {
            "input_otu_id": data.asv_id,
            "in_otu_table": data.asv_id,
            "group_detail": data.group_detail,
            "main_id": str(main_table_id),
            "group": data.group_id,
            "task_id": task_info["task_id"],
            "level": int(9),
            "update_info": json.dumps(update_info)
        }
        if hasattr(data, "method"):
            if data.method in ["none"]:
                options["method"] = ""
            else:
                options["method"] = data.method
        to_file = ["metaasv.export_otu_table_by_detail(in_otu_table)", "metaasv.export_group_table_by_detail(group)"]
        self.set_sheet_data(name=task_name, options=options,
                            main_table_name=main_table_name.strip().split("_")[0] + "/" + main_table_name,
                            params=params_json,
                            module_type=module_type, to_file=to_file)
        task_info = super(Tax4funAction, self).POST()
        task_info["content"] = {"ids": {"id": str(main_table_id), "name": main_table_name}}
        return json.dumps(task_info)

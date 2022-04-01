# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# modified 20180622

import os
import web
import json
import datetime
from bson.objectid import ObjectId
from mainapp.libs.signature import check_sig
from mainapp.models.mongo.dna import Dna
from mainapp.controllers.project.dna_controller import DnaController


class TraitAnnovarAction(DnaController):
    """
    性状分析接口
    """
    def __init__(self):
        super(TraitAnnovarAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        print data
        params = ["trait_id", "task_id", "task_type", "submit_location"]
        for param in params:
            if not hasattr(data, param):
                var = []
                var.append(param)
                info = {"success": False, "info": "缺少%s参数!" % param, "code": "C1700601", "variables": var}
                return json.dumps(info)
        task_result = Dna("dna_gmap").find_one(collection="sg_task", query_dic={"task_id": data.task_id})
        if not task_result:
            var = []
            var.append(data.task_id)
            info = {"success": False, "info": "sg_task表里没有找到任务:%s，请检查!" % data.task_id,
                    "code": "C1700602", "variables": var}
            return json.dumps(info)
        feature_result = Dna("dna_gmap").find_one(collection="sg_feature_file", query_dic={"_id": ObjectId(data.trait_id)})
        if not feature_result:
            var = []
            var.append(data.trait_id)
            info = {"success": False, "info": "sg_feature_file表里没有找到id:%s，请检查!" % data.trait_id,
                    "code": "C1700603", "variables": var}
            return json.dumps(info)
        project_sn = task_result["project_sn"]
        member_id = task_result["member_id"]
        params_json = {
            "trait_id": data.trait_id,
            "task_type": int(data.task_type),
            "submit_location": data.submit_location
        }
        params = json.dumps(params_json, sort_keys=True, separators=(',', ':'))
        main_table_name = 'TraitAnnovar_' + datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        mongo_data = [
            ("name", main_table_name),
            ("status", "start"),
            ("project_sn", project_sn),
            ("task_id", data.task_id),
            ("params", params),
            ("desc", "性状分析"),
            ("created_ts", datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        ]
        main_id = Dna("dna_gmap").insert_main_table(collection="sg_feature", data=mongo_data)
        Dna("dna_gmap").update_db_record(collection="sg_feature", query_dict={"_id": main_id},
                                         update_dict={"main_id": main_id})
        update_info = {str(main_id): "sg_feature"}
        options = {
            "trait_file": data.trait_id,
            "sample_list": data.trait_id,
            "update_info": json.dumps(update_info),
            "main_id": str(main_id),
            "task_id": data.task_id
        }
        to_file = ["dna_gmap.export_trit_file(trait_file)", "dna_gmap.export_sample_file(sample_list)"]
        self.set_sheet_data(name="dna_gmap.trait_annovar", member_id=member_id, project_sn=project_sn, task_id=data.task_id,
                            main_table_name="TraitAnnovar/" + main_table_name, options=options, module_type="tool", params=params,
                            db_type="dna_gmap", analysis_name="TraitAnnovar", to_file=to_file)
        task_info = super(TraitAnnovarAction, self).POST()
        task_info['id'] = str(main_id)
        return json.dumps(task_info)

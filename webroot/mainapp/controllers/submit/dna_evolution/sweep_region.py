# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# modified 20180920

import web
import json
import datetime
from bson.objectid import ObjectId
from mainapp.models.mongo.dna import Dna
from mainapp.libs.signature import check_sig
from mainapp.controllers.project.dna_controller import DnaController


class SweepRegionAction(DnaController):
    """
    受选择区域筛选接口
    """
    def __init__(self):
        super(SweepRegionAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        print data
        params = ["sweep_id", "p_value", "task_id", "task_type", "submit_location", "chongmingming_result"]
        for param in params:
            if not hasattr(data, param):
                info = {"success": False, "info": "缺少%s参数,请检查" % param}
                return json.dumps(info)
        if not hasattr(data, "ld_window_kb") and not hasattr(data, "ld_window_r2"):
            info = {"success": False, "info": "参数%s和%s至少要有一个，请检查" % ("ld_window_kb", "ld_window_r2")}
            return json.dumps(info)
        if float(data.p_value) <= 0 or float(data.p_value) >= 1:
            info = {"success": False, "info": "P值:%s必需为(0,1),请检查" % data.p_value}
            return json.dumps(info)
        if hasattr(data, "ld_window_kb") and float(data.ld_window_kb) <= 0:
            info = {"success": False, "info": "固定大小:%s必须大于0,请检查" % data.ld_window_kb}
            return json.dumps(info)
        if hasattr(data, "ld_window_r2"):
            if float(data.ld_window_r2) <= 0 or float(data.ld_window_r2) >= 1:
                info = {"success": False, "info": "LD选择:%s必须为(0,1),请检查" % data.ld_window_r2}
                return json.dumps(info)
        task_result = Dna("dna_evolution").find_one(collection="sg_task", query_dic={"task_id": data.task_id})
        if not task_result:
            info = {"success": False, "info": "没有在sg_task表里找到task_id:%s,请检查" % data.task_id}
            return json.dumps(info)
        task_id = data.task_id
        project_sn = task_result["project_sn"]
        member_id = task_result["member_id"]
        pop_summary = task_result["pop_summary"]
        sweep_result = Dna("dna_evolution").find_one(collection="sg_sweep", query_dic={"main_id": ObjectId(data.sweep_id)})
        if not sweep_result:
            info = {"success": False, "info": "没有在ssg_sweep表里找到id:%s的结果,请检查" % data.sweep_id}
            return json.dumps(info)
        vcf_file = sweep_result["vcf_path"]
        sweep_dir = sweep_result["sweep_dir"]
        diff_group = ",".join(sweep_result["diff_group"])
        params_json = {
            "sweep_id": data.sweep_id,
            "p_value": data.p_value,
            "submit_location": data.submit_location,
            "chongmingming_result": data.chongmingming_result
        }
        if hasattr(data, "ld_window_kb"):
            params_json["ld_window_kb"] = data.ld_window_kb
        if hasattr(data, "ld_window_r2"):
            params_json["ld_window_r2"] = data.ld_window_r2
        params = json.dumps(params_json, sort_keys=True, separators=(',', ':'))
        # main_table_name = "SweepRegion_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")
        main_table_name = Dna("dna_evolution").set_main_table_name("SweepRegion", data.chongmingming_result)
        mongo_data = [
            ("name", main_table_name),
            ("status", "start"),
            ("project_sn", project_sn),
            ("task_id", task_id),
            ("params", params),
            ("desc", "SweepRegion分析"),
            ("created_ts", datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        ]
        main_id = Dna("dna_evolution").insert_main_table(collection="sg_sweep_region", data=mongo_data)
        Dna("dna_evolution").update_db_record(collection="sg_sweep_region", query_dict={"_id": main_id}, update_dict={"main_id": main_id})
        update_info = {str(main_id): "sg_sweep_region"}
        options = {
            "vcf_file": vcf_file,
            "pop_summary": pop_summary,
            "sweep_dir": sweep_dir,
            "diff_group": diff_group,
            "p_value": data.p_value,
            "update_info": json.dumps(update_info),
            "main_id": str(main_id)
        }
        if hasattr(data, "ld_window_kb"):
            options["ld_window_kb"] = data.ld_window_kb
        if hasattr(data, "ld_window_r2"):
            options["ld_window_r2"] = data.ld_window_r2
        self.set_sheet_data(name="dna_evolution.sweep_region", member_id=member_id, project_sn=project_sn, task_id=task_id,
                            main_table_name="SweepAnalysis/" + main_table_name, options=options, module_type="module", params=params,
                            db_type="dna_evolution", analysis_name="SweepRegion")
        task_info = super(SweepRegionAction, self).POST()
        task_info["id"] = str(main_id)
        return json.dumps(task_info)

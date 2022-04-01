# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# modified 20180913

import web
import json
import datetime
from bson.objectid import ObjectId
from mainapp.models.mongo.dna import Dna
from mainapp.libs.signature import check_sig
from mainapp.controllers.project.dna_controller import DnaController


class SweepAnalysisAction(DnaController):
    """
    选择性消除接口
    """
    def __init__(self):
        super(SweepAnalysisAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        print data
        params = ["vcf_id", "group_dict", "diff_group", "group_id", "analysis_method", "window_size",
                  "task_id", "task_type", "submit_location", "chongmingming_result"]
        for params in params:
            if not hasattr(data, params):
                info = {"success": False, "info": "缺少%s参数" % params}
                return json.dumps(info)
        window_size = float(data.window_size) * 1000
        if window_size <= 0:
            info = {"success": False, "info": "窗口大小:%s必需大于0,请检查" % data.window_size}
            return json.dumps(info)
        if hasattr(data, "max_missing") and data.max_missing:
            if float(data.max_missing) < 0 or float(data.max_missing) > 1:
                info = {"success": False, "info": "缺失率:%s必需为(0,1),请检查" % data.max_missing}
                return json.dumps(info)
        if hasattr(data, "min_maf") and hasattr(data, "max_maf"):
            if data.min_maf and data.max_maf:
                if float(data.min_maf) >= float(data.max_maf):
                    info = {"success": False, "info": "次要等位基因频率前面的值必需小于后面的值:%s-%s,请检查" % (data.min_maf, data.max_maf)}
                    return json.dumps(info)
                if float(data.min_maf) < 0.05:
                    info = {"success": False, "info": "次要等位基因频率:%s必需大于0.05,请检查" % data.min_maf}
                    return json.dumps(info)
                if float(data.max_maf) > 1:
                    info = {"success": False, "info": "次要等位基因频率:%s必需小于1,请检查" % data.max_maf}
                    return json.dumps(info)
        if hasattr(data, "mindp") and hasattr(data, "maxdp"):
            if data.mindp and data.maxdp:
                if int(data.mindp) >= int(data.maxdp):
                    info = {"success": False, "info": "平均测序深度前面的值必需小于后面的值:%s-%s,请检查" % (data.mindp, data.maxdp)}
                    return json.dumps(info)
                if int(data.mindp) < 0:
                    info = {"success": False, "info": "平均测序深度:%s必需大于0,请检查" % data.mindp}
                    return json.dumps(info)
        group_list = json.loads(data.group_dict).keys()
        diff_group = []
        for diff in data.diff_group.split(","):
            group = diff.split("|")[0]
            compare_group = diff.split("|")[1]
            if group not in group_list or compare_group not in group_list:
                info = {"success": False, "info": "差异分组:%s不在分组方案:%s里,请检查" % (diff, group_list)}
                return json.dumps(info)
            diff_group.append(group + "_vs_" + compare_group)
        task_result = Dna("dna_evolution").find_one(collection="sg_task", query_dic={"task_id": data.task_id})
        if not task_result:
            info = {"success": False, "info": "没有在sg_task表里找到task_id:%s,请检查" % data.task_id}
            return json.dumps(info)
        task_id = data.task_id
        project_sn = task_result["project_sn"]
        member_id = task_result["member_id"]
        variant_result = Dna("dna_evolution").find_one(collection="sg_variant_compare_filter", query_dic={"_id": ObjectId(data.vcf_id)})
        if not variant_result:
            info = {"success": False, "info": "没有在sg_variant_compare_filter表里找到id:%s的结果,请检查" % data.vcf_id}
            return json.dumps(info)
        vcf_file = variant_result["vcf_path"]
        # main_table_name = "SweepAnalysis_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")
        main_table_name = Dna("dna_evolution").set_main_table_name("SweepAnalysis", data.chongmingming_result)
        params_json = {
            "vcf_id": data.vcf_id,
            "group_dict": json.loads(data.group_dict),
            "diff_group": data.diff_group,
            "group_id": data.group_id,
            "analysis_method": data.analysis_method,
            "window_size": data.window_size,
            "submit_location": data.submit_location,
            "chongmingming_result": data.chongmingming_result
        }
        if hasattr(data, "max_missing"):
            params_json["max_missing"] = data.max_missing
        if hasattr(data, "mindp"):
            params_json["mindp"] = data.mindp
        if hasattr(data, "maxdp"):
            params_json["maxdp"] = data.maxdp
        if hasattr(data, "min_maf"):
            params_json["min_maf"] = data.min_maf
        if hasattr(data, "max_maf"):
            params_json["max_maf"] = data.max_maf
        params = json.dumps(params_json, sort_keys=True, separators=(',', ':'))
        mongo_data = [
            ("name", main_table_name),
            ("status", "start"),
            ("project_sn", project_sn),
            ("task_id", task_id),
            ("params", params),
            ("diff_group", diff_group),
            ("desc", "SweepAnalysis分析"),
            ("created_ts", datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        ]
        main_id = Dna("dna_evolution").insert_main_table(collection="sg_sweep", data=mongo_data)
        Dna("dna_evolution").update_db_record(collection="sg_sweep", query_dict={"_id": main_id}, update_dict={"main_id": main_id})
        update_info = {str(main_id): "sg_sweep"}
        options = {
            "vcf_file": vcf_file,
            "group_info": data.group_dict,
            "diff_group": ",".join(diff_group),
            "analysis_method": data.analysis_method,
            "window_size": window_size,
            "update_info": json.dumps(update_info),
            "main_id": str(main_id)
        }
        if hasattr(data, "mindp") and data.mindp:
            options["minDP"] = data.mindp
        if hasattr(data, "maxdp") and data.maxdp:
            options["maxDP"] = data.maxdp
        if hasattr(data, "min_maf") and data.min_maf:
            options["min_maf"] = data.min_maf
        if hasattr(data, "max_maf") and data.max_maf:
            options["max_maf"] = data.max_maf
        if hasattr(data, "max_missing") and data.max_missing:
            options["max_missing"] = data.max_missing
        self.set_sheet_data(name="dna_evolution.sweep_analysis", member_id=member_id, project_sn=project_sn, task_id=task_id,
                            main_table_name="SweepAnalysis/" + main_table_name, options=options, module_type="module", params=params,
                            db_type="dna_evolution", analysis_name="SweepAnalysis")
        task_info = super(SweepAnalysisAction, self).POST()
        task_info["id"] = str(main_id)
        return json.dumps(task_info)

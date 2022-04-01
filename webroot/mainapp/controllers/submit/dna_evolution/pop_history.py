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


class PopHistoryAction(DnaController):
    """
    种群历史接口
    """
    def __init__(self):
        super(PopHistoryAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        print data
        params = ["vcf_id", "group_dict", "analysis_method", "max_missing", "chongmingming_result", "group_id",
                  "min_maf", "max_maf", "mindp", "maxdp", "task_id", "task_type", "submit_location"]
        for params in params:
            if not hasattr(data, params):
                info = {"success": False, "info": "缺少%s参数" % params}
                return json.dumps(info)
        if data.analysis_method not in ["psmc", "smc"]:
            info = {"success": False, "info": "分析方法:%s只能是psmc或smc,请检查" % data.analysis_method}
            return json.dumps(info)
        if data.max_missing != "":
            if float(data.max_missing) < 0 or float(data.max_missing) > 1:
                info = {"success": False, "info": "缺失率:%s必需为(0,1),请检查" % data.max_missing}
                return json.dumps(info)
        if data.min_maf != "":
            if float(data.min_maf) < 0.05:
                info = {"success": False, "info": "次要等位基因频率:%s必需大于0.05,请检查" % data.min_maf}
                return json.dumps(info)
        if data.max_maf != "":
            if float(data.max_maf) > 1:
                info = {"success": False, "info": "次要等位基因频率:%s必需小于1,请检查" % data.max_maf}
                return json.dumps(info)
        if data.min_maf != "" and data.max_maf != "":
            if float(data.min_maf) >= float(data.max_maf):
                info = {"success": False, "info": "次要等位基因频率前面的值必需小于后面的值:%s-%s,请检查" % (data.min_maf, data.max_maf)}
                return json.dumps(info)
        if data.mindp != "":
            if int(data.mindp) < 0:
                info = {"success": False, "info": "平均测序深度:%s必需大于0,请检查" % data.mindp}
                return json.dumps(info)
            if data.maxdp != "":
                if int(data.mindp) >= int(data.maxdp):
                    info = {"success": False, "info": "平均测序深度前面的值必需小于后面的值:%s-%s,请检查" % (data.mindp, data.maxdp)}
                    return json.dumps(info)
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
        main_table_name = Dna("dna_evolution").set_main_table_name("PopHistory", data.chongmingming_result)
        params_json = {
            "vcf_id": data.vcf_id,
            "group_dict": json.loads(data.group_dict),
            "group_id": data.group_id,
            "analysis_method": data.analysis_method,
            "max_missing": data.max_missing,
            "min_maf": data.min_maf,
            "max_maf": data.max_maf,
            "mindp": data.mindp,
            "maxdp": data.maxdp,
            "submit_location": data.submit_location,
            "chongmingming_result": data.chongmingming_result
        }
        params = json.dumps(params_json, sort_keys=True, separators=(',', ':'))
        # main_table_name = "PopHistory_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")
        mongo_data = [
            ("name", main_table_name),
            ("status", "start"),
            ("project_sn", project_sn),
            ("task_id", task_id),
            ("params", params),
            ("desc", "PopHistory分析"),
            ("created_ts", datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        ]
        main_id = Dna("dna_evolution").insert_main_table(collection="sg_history", data=mongo_data)
        Dna("dna_evolution").update_db_record(collection="sg_history", query_dict={"_id": main_id}, update_dict={"main_id": main_id})
        update_info = {str(main_id): "sg_history"}
        options = {
            "vcf_file": vcf_file,
            "group_info": data.group_dict,
            "analysis_method": data.analysis_method,
            "update_info": json.dumps(update_info),
            "main_id": str(main_id)
        }
        if data.max_missing != "":
            options["max_missing"] = data.max_missing
        if data.min_maf != "":
            options["min_maf"] = data.min_maf
        if data.max_maf != "":
            options["max_maf"] = data.max_maf
        if data.mindp != "":
            options["minDP"] = data.mindp
        if data.maxdp != "":
            options["maxDP"] = data.maxdp
        self.set_sheet_data(name="dna_evolution.pop_history", member_id=member_id, project_sn=project_sn, task_id=task_id,
                            main_table_name="PopHistory/" + main_table_name, options=options, module_type="module", params=params,
                            db_type="dna_evolution", analysis_name="PopHistory")
        task_info = super(PopHistoryAction, self).POST()
        task_info["id"] = str(main_id)
        return json.dumps(task_info)

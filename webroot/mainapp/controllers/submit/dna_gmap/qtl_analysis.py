# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# modified 20180705

import os
import web
import json
import datetime
from bson.objectid import ObjectId
from mainapp.libs.signature import check_sig
from mainapp.models.mongo.dna import Dna
from mainapp.controllers.project.dna_controller import DnaController


class QtlAnalysisAction(DnaController):
    """
    QTL定位分析接口
    """
    def __init__(self):
        super(QtlAnalysisAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        print data
        params = ["lg_id", "trait_id", "location_method", "pm_num", "task_id", "task_type", "submit_location"]
        for param in params:
            if not hasattr(data, param):
                var = []
                var.append(param)
                info = {"success": False, "info": "缺少%s参数!" % param, "code": "C1700401", "variables": var}
                return json.dumps(info)
        if not data.p_value and not data.lod_value:
            info = {"success": False, "info": "缺少阈值参数!", "code": "C1700402", "variables": ""}
            return json.dumps(info)
        if data.location_method not in ["cim", "scanone"]:
            var = []
            var.append(data.location_method)
            info = {"success": False, "info": "定位方法只能为cim/scanone,而不能是:%s" % data.location_method,
                    "code": "C1700403", "variables": var}
            return json.dumps(info)
        if data.pm_num < 0:
            info = {"success": False, "info": "置换检验次数只能是大于等于0的正整数", "code": "C1700404", "variables": ""}
            return json.dumps(info)
        task_result = Dna("dna_gmap").find_one(collection="sg_task", query_dic={"task_id": data.task_id})
        if not task_result:
            var = []
            var.append(data.task_id)
            info = {"success": False, "info": "sg_task表里没有找到任务:%s，请检查!" % data.task_id,
                    "code": "C1700405", "variables": var}
            return json.dumps(info)
        project_sn = task_result["project_sn"]
        member_id = task_result["member_id"]
        pop_type = "F1" if task_result["poptype"] == "CP" else task_result["poptype"]
        pid = task_result["pid"]
        mid = task_result["mid"]
        # pop_final_vcf = self.get_path(task_result["pop_final_vcf"])
        # pop_summary = self.get_path(task_result["pop_summary"])
        pop_final_vcf = Dna("dna_gmap").set_file_path(data.task_id, task_result["pop_final_vcf"], data.client)
        pop_summary = Dna("dna_gmap").set_file_path(data.task_id, task_result["pop_summary"], data.client)
        if pop_type == "F1" and data.location_method == "cim":
            info = {"success": False, "info": "群体类型为F1的时候定位方法不能是cim,请检查", "code": "C1700406", "variables": ""}
            return json.dumps(info)
        lg_result = Dna("dna_gmap").find_one(collection="sg_lg", query_dic={"main_id": ObjectId(data.lg_id)})
        if not lg_result:
            var = []
            var.append(data.lg_id)
            info = {"success": False, "info": "sg_lg里没有找到id:%s对应的结果，请检查" % data.lg_id,
                    "code": "C1700407", "variables": var}
            return json.dumps(info)
        # if pop_type == "F1":
        #     # marker_info_path = self.get_path(lg_result["marker_info_path"])
        #     marker_info_path = Dna("dna_gmap").set_file_path(data.task_id, lg_result["marker_info_path"], data.client)
        # else:
        #     # total_map = self.get_path(lg_result["total_csv_path"])
        #     # total_phase = self.get_path(lg_result["total_loc_path"])
        #     total_map = Dna("dna_gmap").set_file_path(data.task_id, lg_result["total_csv_path"], data.client)
        #     total_phase = Dna("dna_gmap").set_file_path(data.task_id, lg_result["total_loc_path"], data.client)
        if pop_type != "F1":
            total_phase = Dna("dna_gmap").set_file_path(data.task_id, lg_result["total_loc_path"], data.client)
        marker_info_path = Dna("dna_gmap").set_file_path(data.task_id, lg_result["marker_info_path"], data.client)
        feature_result = Dna("dna_gmap").find_one(collection="sg_feature_file", query_dic={"_id": ObjectId(data.trait_id)})
        if not feature_result:
            var = []
            var.append(data.trait_id)
            info = {"success": False, "info": "sg_feature_file里没有找到id:%s对应的结果，请检查" % data.trait_id,
                    "code": "C1700408", "variables": var}
            return json.dumps(info)
        params_json = {
            "lg_id": data.lg_id,
            "trait_id": data.trait_id,
            "location_method": data.location_method,
            "pm_num": data.pm_num,
            "task_type": int(data.task_type),
            "submit_location": data.submit_location
        }
        if data.p_value:
            params_json["p_value"] = data.p_value
        else:
            params_json["lod_value"] = data.lod_value
        params = json.dumps(params_json, sort_keys=True, separators=(',', ':'))
        tm = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        main_table_name = "QtlAnalysis_" + tm
        mongo_data = [
            ("name", main_table_name),
            ("status", "start"),
            ("project_sn", project_sn),
            ("task_id", data.task_id),
            ("params", params),
            ("desc", "QTL定位分析"),
            ("created_ts", tm)
        ]
        main_id = Dna("dna_gmap").insert_main_table(collection="sg_qtl", data=mongo_data)
        Dna("dna_gmap").update_db_record(collection="sg_qtl", query_dict={"_id": main_id},
                                         update_dict={"main_id": main_id})
        update_info = {str(main_id): "sg_qtl"}
        options = {
            "trait_dir": data.trait_id,
            "pop_final_vcf": pop_final_vcf,
            "pop_summary": pop_summary,
            "pop_type": pop_type,
            "location_method": data.location_method,
            "pm_num": data.pm_num,
            "pid": pid,
            "mid": mid,
            "update_info": json.dumps(update_info),
            "main_id": str(main_id)
        }
        if data.p_value:
            options["p_value"] = data.p_value
        else:
            options["lod_value"] = data.lod_value
        # if pop_type == "F1":
        #     options["total_phase"] = marker_info_path
        #     options["total_map"] = marker_info_path
        # else:
        #     options["total_map"] = total_map
        if pop_type == "F1":
            options["total_phase"] = marker_info_path
        options["total_map"] = marker_info_path
        to_file = ["dna_gmap.export_trit_dir(trait_dir)"]
        self.set_sheet_data(name="dna_gmap.qtl_analysis", member_id=member_id, project_sn=project_sn, task_id=data.task_id,
                            main_table_name="QtlAnalysis/" + main_table_name, options=options, module_type="module", params=params,
                            db_type="dna_gmap", analysis_name="QtlAnalysis", to_file=to_file)
        task_info = super(QtlAnalysisAction, self).POST()
        task_info['id'] = str(main_id)
        return json.dumps(task_info)

    def get_path(self, rere_path):
        if os.path.exists(rere_path):
            return rere_path
        target_path = os.path.join("s3://", rere_path)
        return target_path

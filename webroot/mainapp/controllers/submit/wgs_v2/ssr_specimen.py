# -*- coding: utf-8 -*-
# __author__: zengjing
# modified: 20190321

import web
import json
import datetime
from mainapp.models.mongo.dna import Dna
from mainapp.libs.signature import check_sig
from mainapp.controllers.project.dna_controller import DnaController


class SsrSpecimenAction(DnaController):
    """
    SSR标记开发接口
    """
    def __init__(self):
        super(SsrSpecimenAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        print data
        params = ["task_id", "task_type", "submit_location", "chongmingming_result"]
        for param in params:
            if not hasattr(data, param):
                info = {"success": False, "info": "缺少参数: %s" % param}
                return json.dumps(info)
        result = Dna("dna_wgs_v2").find_one(collection="sg_task", query_dic={"task_id": data.task_id})
        try:
            project_sn = result["project_sn"]
            member_id = result["member_id"]
            bam_path = result["bam_path"]
            genome_version_id = result["genome_version_id"]
        except:
            info = {"success": False, "info": "sg_task表没有找到信息,请检查"}
            return json.dumps(info)
        results = Dna("dna_wgs_v2").find_many(collection="sg_specimen", query_dic={"task_id": data.task_id})
        try:
            specimen_list = []
            for result in results:
                specimen_list.append(result["analysis_name"])
            specimen_list = list(set(specimen_list))
        except:
            info = {"success": False, "info": "获取样本列表失败,请检查"}
            return json.dumps(info)
        result = Dna("dna_wgs_v2").find_one(collection="sg_species_version", query_dic={"_id": genome_version_id})
        try:
            ref_ssr = result["ssr_path"]
            total_chrlist = result["total_chrlist"]
        except:
            info = {"success": False, "info": "sg_species_version表没有找到ssr信息信息,请检查"}
            return json.dumps(info)
        # if data.client == "client03":
        #     ref_misa = "/mnt/ilustre/users/sanger-dev/app/database/dna_geneome/" + ref_ssr + "ssr.ref.result.xls"
        #     ref_chrlist = "/mnt/ilustre/users/sanger-dev/app/database/dna_geneome/" + total_chrlist
        # else:
        #     ref_misa = "/mnt/lustre/users/sanger/app/database/dna_geneome/" + ref_ssr + "ssr.ref.result.xls"
        #     ref_chrlist = "/mnt/lustre/users/sanger/app/database/dna_geneome/" + total_chrlist
        ref_misa = ref_ssr + "ssr.ref.result.xls"
        ref_chrlist = total_chrlist
        params_json = {
            "task_type": int(data.task_type),
            "task_id": data.task_id,
            "submit_location": data.submit_location,
            "chongmingming_result": data.chongmingming_result
        }
        params = json.dumps(params_json, sort_keys=True, separators=(",", ":"))
        if hasattr(data, "chongmingming_result") and data.chongmingming_result:
            main_table_name = data.chongmingming_result
        else:
            main_table_name = "SSR_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:3]
        mongo_data = [
            ("project_sn", project_sn),
            ("task_id", data.task_id),
            ("params", params),
            ("specimen_list", specimen_list),
            ("name", main_table_name),
            ("desc", "SSR标记分析"),
            ("status", "start"),
            ("created_ts", datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        ]
        main_id = Dna("dna_wgs_v2").insert_main_table(collection="sg_ssr_marker", data=mongo_data)
        Dna("dna_wgs_v2").update_db_record(collection="sg_ssr_marker", query_dict={"_id": main_id}, update_dict={"main_id": main_id})
        update_info = {str(main_id): "sg_ssr_marker"}
        options = {
            "bam_path": bam_path,
            "ref_chrlist": ref_chrlist,
            "ref_misa": ref_misa,
            "update_info": json.dumps(update_info),
            "main_id": str(main_id),
            "task_id": data.task_id
        }
        main_table_name = "SSR/SSR_" + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3])
        self.set_sheet_data(name="wgs_v2.report.ssr_specimen", member_id=member_id, project_sn=project_sn, task_id=data.task_id,
                            main_table_name=main_table_name, options=options, module_type="workflow", params=params, db_type="dna_wgs_v2")
        task_info = super(SsrSpecimenAction, self).POST()
        task_info["id"] = str(main_id)
        return json.dumps(task_info)

# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# modified 2018.04.24

import os
import web
import json
import datetime
from bson.objectid import ObjectId
from mainapp.libs.signature import check_sig
from mainapp.models.mongo.dna import Dna
from mainapp.controllers.project.dna_controller import DnaController


class PrimerDesignAction(DnaController):
    """
    引物设计接口
    """
    def __init__(self):
        super(PrimerDesignAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        print data
        params = ["set_id", "tm1", "tm2", "product_size", "primer_num", "task_id", "task_type", "submit_location"]
        for param in params:
            if not hasattr(data, param):
                var = []
                var.append(param)
                info = {"success": False, "info": "缺少%s参数!" % param, "code": "C3100801", "variables": var}
                return json.dumps(info)
        set_id = data.set_id
        if float(data.tm1) >= float(data.tm2):
            var = []
            var.append(data.tm1)
            var.append(data.tm2)
            info = {"success": False, "info": "tm1:%s不能大于等于tm2:%s,请检查!" % (data.tm1, data.tm2),
                    "code": "C3100802", "variables": var}
            return json.dumps(info)
        if int(data.primer_num) > 5 or int(data.primer_num) < 1:
            var = []
            var.append(data.primer_num)
            info = {"success": False, "info": "primer_num:%s范围只能是[2,5],请检查!" % data.primer_num,
                    "code": "C3100803", "variables": var}
            return json.dumps(info)
        if not hasattr(data, "project_type"):
            db = "dna_wgs"
        elif data.project_type == "dna_gmap":
            db = "dna_gmap"
        else:
            db = data.project_type
        table_name = ""
        if set_id != "all":
            result = Dna(db).find_one(collection="sg_site_set", query_dic={"task_id": data.task_id, "set_id": set_id})
            if not result:
                var = []
                var.append(data.primer_num)
                info = {"success": False, "info": "sg_site_set表里没有该set_id: %s,请检查!" % set_id,
                        "code": "C3100804", "variables": var}
                return json.dumps(info)
            type = result["type"]
            table_name = "_".join(result["name"].strip().split("_")[:-2])
            if type.startswith("snp"):
                compare_result = Dna(db).find_one(collection="sg_snp_compare", query_dic={"main_id": ObjectId(set_id)})
            else:
                compare_result = Dna(db).find_one(collection="sg_indel_compare", query_dic={"main_id": ObjectId(set_id)})
            if not compare_result:
                var = []
                var.append(set_id)
                info = {"success": False, "info": "sg_snp_compare/sg_indel_compare表里没有该id:%s,请检查!" % set_id,
                        "code": "C3100805", "variables": var }
                return json.dumps(info)
            if "download_path" in compare_result.keys():
                # diff_variant = os.path.join("s3://", compare_result["download_path"])
                diff_variant = Dna(db).set_file_path(data.task_id, compare_result["download_path"], data.client)
            else:
                diff_variant = compare_result["diff_variant"]
        else:
            table_name = "PrimerDesign"
            task_result = Dna(db).find_one(collection="sg_task", query_dic={"task_id": data.task_id})
            if not task_result:
                var = []
                var.append(set_id)
                info = {"success": False, "info": "sg_task表里没有该task_id:%s,请检查!" % data.task_id,
                        "code": "C3100806", "variables": var}
                return json.dumps(info)
            diff_variant = task_result["pop_final_vcf"]
        diff_variant = Dna(db).set_file_path(data.task_id, diff_variant, data.client)
        task_result = Dna(db).find_one(collection="sg_task", query_dic={"task_id": data.task_id})
        if not task_result:
            var = []
            var.append(data.task_id)
            info = {"success": False, "info": "sg_task表里没有该task_id:%s,请检查!" % data.task_id,
                    "code": "C3100807", "variables": var}
            return json.dumps(info)
        project_sn = task_result["project_sn"]
        member_id = task_result["member_id"]
        genome_version_id = task_result["genome_version_id"]
        ref_result = Dna(db).find_one(collection="sg_species_version", query_dic={"_id": genome_version_id})
        # noinspection PyBroadException
        try:
            ref_fa = ref_result['ref']
        except:
            info = {"success": False, "info": "sg_species_version表里没有ref信息，请检查!",
                    "code": "C3100808", "variables": ""}
            return json.dumps(info)
        params_json = {
            "set_id": data.set_id,
            "tm1": data.tm1,
            "tm2": data.tm2,
            "product_size": data.product_size,
            "primer_num": data.primer_num,
            "task_type": data.task_type,
            "submit_location": data.submit_location
        }
        params = json.dumps(params_json, sort_keys=True, separators=(',', ':'))
        mongo_data = [
            ("name", "PrimerDesign_" + table_name + "_" + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S"))),
            ("status", "start"),
            ("project_sn", project_sn),
            ("task_id", data.task_id),
            ("params", params),
            ("desc", "引物设计分析主表"),
            ("created_ts", datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        ]
        main_id = Dna(db).insert_main_table(collection="sg_primer", data=mongo_data)
        Dna(db).update_db_record(collection="sg_primer", query_dict={"_id": main_id},
                                 update_dict={"main_id": main_id})
        update_info = {str(main_id): "sg_primer"}
        options = {
            "ref_fa": ref_fa,
            "diff_variant": diff_variant,
            "tm1": float(data.tm1),
            "tm2": float(data.tm2),
            "product_size": data.product_size,
            "primer_num": int(data.primer_num),
            "update_info": json.dumps(update_info),
            "main_id": str(main_id)
        }
        if hasattr(data, "project_type"):
            options["project_type"] = data.project_type
        main_table_name = "PrimerDesign/PrimerDesign_" + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3])
        self.set_sheet_data(name="wgs.report.primer_design", member_id=member_id, project_sn=project_sn, task_id=data.task_id,
                            main_table_name=main_table_name, options=options, module_type="workflow", params=params, db_type=db)
        task_info = super(PrimerDesignAction, self).POST()
        task_info['id'] = str(main_id)
        return json.dumps(task_info)

# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# modified 2018.05.18

import os
import web
import json
import datetime
from bson.objectid import ObjectId
from mainapp.libs.signature import check_sig
from mainapp.models.mongo.dna import Dna
from mainapp.controllers.project.dna_controller import DnaController


class SsrCompareAction(DnaController):
    """
    SSR比较分析接口
    """
    def __init__(self):
        super(SsrCompareAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        print data
        params = ["ssr_specimen_id1", "ssr_specimen_id2", "sample1", "sample2", "is_same", "task_id", "task_type", "submit_location"]
        for param in params:
            if not hasattr(data, param):
                var = []
                var.append(param)
                info = {"success": False, "info": "缺少%s参数!" % param, "code": "C3101201", "variables": var}
                return json.dumps(info)
        if data.is_same not in ["true", "false"]:
            var = []
            var.append(data.is_same)
            info = {"success": False, "info": "相同与不同的类型%s不合法!, 必须为true或者false" % data.is_same,
                    "code": "C3101202", "variables": var}
            return json.dumps(info)
        if not hasattr(data, "project_type"):
            db = "dna_wgs"
        elif data.project_type == "dna_gmap":
            db = "dna_gmap"
        else:
            db = data.project_type
        result = Dna(db).find_one(collection="sg_task", query_dic={"task_id": data.task_id})
        try:
            project_sn = result["project_sn"]
            member_id = result["member_id"]
            genome_version_id = result["genome_version_id"]
        except:
            info = {"success": False, "info": "sg_task表里没有project_sn, member_id信息，请检查!",
                    "code": "C3101203", "variables": ""}
            return json.dumps(info)
        ref_result = Dna(db).find_one(collection="sg_species_version", query_dic={"_id": genome_version_id})
        try:
            ssr_path = ref_result["ssr_path"]
            if data.client == "client03":
                ref_db = os.path.join("/mnt/ilustre/users/sanger-dev/app/database/dna_wgs_geneome/", ssr_path + "makedbblast/ref")
            else:
                ref_db = os.path.join("/mnt/lustre/users/sanger/app/database/dna_wgs_geneome/", ssr_path + "makedbblast/ref")
        except:
            info = {"success": False, "info": "sg_species_version表里没有ref、gff信息，请检查!",
                    "code": "C3101204", "variables": ""}
            return json.dumps(info)
        ssr_result_1 = Dna(db).find_one(collection="sg_ssr_specimen", query_dic={"main_id": ObjectId(data.ssr_specimen_id1)})
        try:
            # output_dir = os.path.join("s3://", ssr_result_1["download_path"])
            output_dir = Dna(db).set_file_path(data.task_id, ssr_result_1["download_path"], data.client)
            query_fa1 = os.path.join(output_dir, data.sample1 + ".denovo.scafSeq")
            ssr_result1 = os.path.join(output_dir, data.sample1 + ".result")
            if not os.path.exists(query_fa1):
                var = []
                var.append(data.sample1)
                info = {"success": False, "info": "样本%s的scafSeq文件不存在，请检查!" % (data.sample1),
                        "code": "C3101205", "variables": var}
                return json.dumps(info)
            if not os.path.exists(ssr_result1):
                var = []
                var.append(data.sample1)
                info = {"success": False, "info": "样本%s的result文件不存在，请检查!" % (data.sample1),
                        "code": "C3101206", "variables": var}
                return json.dumps(info)
        except:
            var = []
            var.append(data.ssr_specimen_id1)
            info = {"success": False, "info": "sg_ssr_specimen没找到%s对应的output_dir信息，请检查!" % (data.ssr_specimen_id1),
                    "code": "C3101207", "variables": var}
            return json.dumps(info)
        ssr_result_2 = Dna(db).find_one(collection="sg_ssr_specimen", query_dic={"main_id": ObjectId(data.ssr_specimen_id2)})
        try:
            # output_dir = os.path.join("s3://", ssr_result_2["download_path"])
            output_dir = Dna(db).set_file_path(data.task_id, ssr_result_2["download_path"], data.client)
            query_fa2 = os.path.join(output_dir, data.sample2 + ".denovo.scafSeq")
            ssr_result2 = os.path.join(output_dir, data.sample2 + ".result")
            if not os.path.exists(query_fa2):
                var = []
                var.append(data.sample2)
                info = {"success": False, "info": "样本%s的scafSeq文件不存在，请检查!" % (data.sample2),
                        "code": "C3101208", "variables": var}
                return json.dumps(info)
            if not os.path.exists(ssr_result2):
                var = []
                var.append(data.sample2)
                info = {"success": False, "info": "样本%s的result文件不存在，请检查!" % (data.sample2),
                        "code": "C3101209", "variables": var}
                return json.dumps(info)
        except:
            var = []
            var.append(data.ssr_specimen_id2)
            info = {"success": False, "info": "sg_ssr_specimen没找到%s对应的output_dir信息，请检查!" % (data.ssr_specimen_id2),
                    "code": "C3101210", "variables": var}
            return json.dumps(info)
        params_json = {
            "ssr_specimen_id1": data.ssr_specimen_id1,
            "ssr_specimen_id2": data.ssr_specimen_id2,
            "sample1": data.sample1,
            "sample2": data.sample2,
            "is_same": data.is_same,
            "task_type": int(data.task_type),
            "submit_location": data.submit_location
        }
        params = json.dumps(params_json, sort_keys=True, separators=(',', ':'))
        main_table_name = 'SsrCompare_' + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        mongo_data = [
            ("name", main_table_name),
            ("status", "start"),
            ("project_sn", project_sn),
            ("task_id", data.task_id),
            ("params", params),
            ("desc", "SSR比较分析！"),
            ("created_ts", datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        ]
        main_id = Dna(db).insert_main_table(collection="sg_ssr_compare", data=mongo_data)
        Dna(db).update_db_record(collection="sg_ssr_compare", query_dict={"_id": main_id},
                                 update_dict={"main_id": main_id})
        update_info = {str(main_id): "sg_ssr_compare"}
        options = {
            "query_fa1": query_fa1,
            "query_fa2": query_fa2,
            "ssr_result1": ssr_result1,
            "ssr_result2": ssr_result2,
            "dbname_nsq": ref_db,
            "is_same": data.is_same,
            "update_info": json.dumps(update_info),
            "main_id": str(main_id)
        }
        if hasattr(data, "project_type"):
            options["project_type"] = data.project_type
        main_table_name = "SsrCompare/SsrCompare_" + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3])
        self.set_sheet_data(name="wgs.report.ssr_compare", member_id=member_id, project_sn=project_sn, task_id=data.task_id,
                            main_table_name=main_table_name, options=options, module_type="workflow", params=params, db_type=db)
        task_info = super(SsrCompareAction, self).POST()
        task_info['id'] = str(main_id)
        return json.dumps(task_info)

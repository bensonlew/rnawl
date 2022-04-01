# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# modified 20180621

import os
import web
import json
import datetime
from bson.objectid import ObjectId
from mainapp.libs.signature import check_sig
from mainapp.models.mongo.dna import Dna
from mainapp.controllers.project.dna_controller import DnaController


class BinnerCalculateAction(DnaController):
    """
    Bin marker分析接口
    """
    def __init__(self):
        super(BinnerCalculateAction, self).__init__(instant=False)
    @check_sig
    def POST(self):
        data = web.input()
        print data
        params = ["matrix_id", "window_size", "window_step", "upload_marker",
                  "task_id", "task_type", "submit_location"]
        for param in params:
            if not hasattr(data, param):
                var = []
                var.append(param)
                info = {"success": False, "info": "缺少%s参数!" % param, "code": "C1700101", "variables": var}
                return json.dumps(info)
        task_result = Dna("dna_gmap").find_one(collection="sg_task", query_dic={"task_id": data.task_id})
        if not task_result:
            var = []
            var.append(data.task_id)
            info = {"success": False, "info": "sg_task表中没有{}对应的信息!" % data.task_id, "code": "C1700102", "variables": var}
            return json.dumps(info)
        else:
            # noinspection PyBroadException
            try:
                project_sn = task_result["project_sn"]
                member_id = task_result["member_id"]
                pop_final_vcf = task_result["pop_final_vcf"]
                pop_type = task_result["poptype"]
            except:
                info = {"success": False, "info": "sg_task表中没有project_sn, member_id, pop_final_vcf, poptype!", "code" : "C1700110"}
                return json.dumps(info)
        pop_type = "F1" if pop_type == "CP" else pop_type
        pop_final_vcf = Dna("dna_gmap").set_file_path(data.task_id, pop_final_vcf, data.client)
        # pop_final_vcf = self.checkout_file(task_result["pop_final_vcf"], data.client)
        result = Dna("dna_gmap").find_one(collection="sg_marker", query_dic={"main_id": ObjectId(data.matrix_id)})
        if not result:
            var = []
            var.append(data.matrix_id)
            info = {"success": False, "info": "sg_marker表里没有matrix_id:%s对应的信息，请检查!" % data.matrix_id,
                    "code": "C1700103", "variables": var}
            return json.dumps(info)
        else:
            # noinspection PyBroadException
            try:
                genotype_matrix = result["filtered_marker_path"]
            except:
                info = {"success": False, "info": "sg_marker表中没有filtered_marker_path!", "code" : "C1700111"}
                return json.dumps(info)
        genotype_matrix = Dna("dna_gmap").set_file_path(data.task_id, genotype_matrix, data.client)
        params = json.loads(result["params"])
        if "upload_marker" in params:
            upload_marker = params["upload_marker"]
        else:
            if data.upload_marker == "1":

                info = {"success": False, "info": "生成此分型矩阵的时候没有上传标记列表，因此'是否选择标记列表做Bin'不能选择'是'!",
                        "code": "C1700104", "variables": ""}
                return json.dumps(info)
        params_json = {
            "matrix_id": data.matrix_id,
            "window_size": data.window_size,
            "window_step": data.window_step,
            "upload_marker": data.upload_marker,
            "task_type": int(data.task_type),
            "submit_location": data.submit_location
        }
        params = json.dumps(params_json, sort_keys=True, separators=(',', ':'))
        main_table_name = 'BinMarker_' + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")
        mongo_data = [
            ("name", main_table_name),
            ("status", "start"),
            ("project_sn", project_sn),
            ("task_id", data.task_id),
            ("params", params),
            ("desc", "Bin Marker分析"),
            ("created_ts", datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        ]
        main_id = Dna("dna_gmap").insert_main_table(collection="sg_binmarker", data=mongo_data)
        Dna("dna_gmap").update_db_record(collection="sg_binmarker", query_dict={"_id": main_id},
                                         update_dict={"main_id": main_id})
        update_info = {str(main_id): "sg_binmarker"}
        options = {
            "pop_final_vcf": pop_final_vcf,
            "genotype_matrix": genotype_matrix,
            "pop_type": pop_type,
            "window_size": data.window_size,
            "window_step": data.window_step,
            "update_info": json.dumps(update_info),
            "main_id": str(main_id)
        }
        if data.upload_marker == "1":
            options["marker_upload"] = upload_marker
        self.set_sheet_data(name="dna_gmap.report.binner_calculate", member_id=member_id, project_sn=project_sn, task_id=data.task_id,
                            main_table_name="BinMarker/" + main_table_name, options=options, module_type="workflow", params=params,
                            db_type="dna_gmap", analysis_name="BinMarker")
        task_info = super(BinnerCalculateAction, self).POST()
        task_info['id'] = str(main_id)
        return json.dumps(task_info)

    def checkout_file(self, file_path, client):
        if os.path.exists(file_path):
            path = file_path
        elif client == 'client01':
            path = os.path.join("/mnt/ilustre/data", file_path)
        else:
            path = os.path.join("/mnt/ilustre/tsanger-data", file_path)
        if not os.path.exists(path):
            var = []
            var.append(path)
            info = {"success": False, "info": "文件:%s不存在，请检查!" % path, "code": "C1700105", "variables": var}
            return json.dumps(info)
        return path

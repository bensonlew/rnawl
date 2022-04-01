# -*- coding: utf-8 -*-
# __author__ = 'qing_mei
# modified 20180626
# modified 20180707
# controller.submit

import web
import json
import datetime
import os
from bson.objectid import ObjectId
from types import StringTypes
from mainapp.libs.signature import check_sig
from mainapp.models.mongo.dna import Dna
from mainapp.controllers.project.dna_controller import DnaController


class MarkerFilterAction(DnaController):
    """
    遗传图谱 遗传标记筛选分析接口
    params通过导表的
    """

    def __init__(self):
        super(MarkerFilterAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()  # 必传client
        print data
        params = ["type", "pdep", "odep", "miss_tatio", "signif", "matrix_id", "marker_upload", "generationid",
                  "submit_location"]
        for param in params:
            if not hasattr(data, param):
                var = []
                var.append(param)
                info = {"success": False, "info": "缺少%s参数!" % param, "code": "C1700301", "variables": var}
                return json.dumps(info)
        if data.matrix_id == "":
            var = []
            var.append(data.matrix_id)
            info = {"success": False, "info": "%s参数为空，请核查!" % data.matrix_id,
                    "code": "C1700302", "variables": var}
            return json.dumps(info)
        task_result = Dna("dna_gmap").find_one(collection="sg_task", query_dic={"task_id": data.task_id})
        if not task_result:
            var = []
            var.append(data.task_id)
            info = {"success": False,
                    "info": "sg_task表里没有task_id: %s对应的信息，请检查!" % data.task_id,
                    "code": "C1700303", "variables": var}
            return json.dumps(info)
        project_sn = task_result["project_sn"]
        member_id = task_result["member_id"]
        popt = task_result['poptype']
        option_popt = task_result['poptype']    # 传递workflow用
        if popt.lower() == 'ril':   # 表内为RiL
            popt = popt + str(task_result['daishu'])
            option_popt = 'ri' + str(task_result['daishu'])
        params_json = {
            "pid": task_result['pid'],
            "mid": task_result['mid'],
            "type": data.type if data.type else "ALL",
            "pdep": data.pdep,
            "odep": data.odep,
            "popt": popt,
            "miss_tatio": float(data.miss_tatio),
            "signif": float(data.signif),
            "matrix_id": data.matrix_id,    # 页面传分型矩阵sg_subtype_matrix的_id—20180724
            "submit_location": data.submit_location
        }
        if data.marker_upload != "":  # 老师上传文件
            params_json['marker_upload'] = data.marker_upload.lstrip('/')  # 这边周勇那边拼接好具体的路径传过来
            # params_json['marker_upload'] = "s3://{}".format(data.marker_upload.lstrip('/'))
        if data.generationid:
            params_json['generationid'] = data.generationid
        # base_path = 's3://'
        matrix_result = Dna("dna_gmap").find_one(collection="sg_subtype_matrix", query_dic={"_id": ObjectId(data.matrix_id)})
        detail_info_path = matrix_result['path']
        # detail_info_path = os.path.join(base_path, detail_info_path)
        detail_info_path = Dna("dna_gmap").set_file_path(data.task_id, detail_info_path, data.client)
        if matrix_result['origin_id'] != '':
            sg_marker_id = matrix_result['origin_id']
            master_result = Dna("dna_gmap").find_one(collection="sg_marker", query_dic={"_id": ObjectId(sg_marker_id)})
            master_params_json = json.loads(master_result['params'])
            for i in ['pdep', 'odep', 'miss_tatio', 'signif']:
                if params_json[i] == "":
                    params_json[i] = master_params_json[i]
        if params_json['miss_tatio'] == "":     # 书写原因,pl程序默认数值。我需要传入主表内，因此在coltroller里写
            params_json['miss_tatio'] = 30
        if params_json['signif'] == "":
            params_json['signif'] = 0.05    # 同上
        params = json.dumps(params_json, sort_keys=True, separators=(',', ':'))     # 返回json
        main_table_name = 'MarkerFilter_' + \
            datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")
        mongo_data = [
            ("name", main_table_name),
            ("project_sn", project_sn),
            ("task_id", data.task_id),
            ("created_ts", datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ("type", data.type),
            # ("chr_list", chr_list),
            ("status", "start"),
            ("params", params),
            ("desc", "遗传标记筛选"),
        ]
        main_id = Dna("dna_gmap").insert_main_table(collection="sg_marker", data=mongo_data)
        Dna("dna_gmap").update_db_record(collection="sg_marker", query_dict={"_id": main_id},
                                         update_dict={"main_id": main_id})
        update_info = {str(main_id): "sg_marker"}
        result = Dna("dna_gmap").find_one(collection="sg_marker", query_dic={"_id": main_id})
        if not result:
            var = []
            var.append(main_id)
            info = {"success": False,
                    "info": "sg_marker表里没有main_id: %s对应的信息，请检查!" % main_id,
                    "code": "C1200404", "variables": var}
            return json.dumps(info)
        types = ''
        if len(data.type.strip(',').split(',')) != 2:
            types = data.type.strip(',').upper()
        else:
            types = 'ALL'
        # vcf_path = os.path.join(base_path, task_result['pop_final_vcf'])
        vcf_path = Dna("dna_gmap").set_file_path(data.task_id, task_result['pop_final_vcf'], data.client)
        # ref_chrlist_path = os.path.join(base_path, task_result['ref_chrlist'])
        ref_chrlist_path = Dna("dna_gmap").set_file_path(data.task_id, task_result['ref_chrlist'], data.client)
        options = {
            "main_id": str(main_id),
            "task_id": data.task_id,
            "vcf": vcf_path,
            "detail_info": detail_info_path,
            "type": types,
            "pdep": data.pdep,
            "odep": data.odep,
            "popt": option_popt,
            "miss_tatio": float(data.miss_tatio) / 100,
            "signif": float(data.signif),
            "update_info": json.dumps(update_info),     # 给work_flow用
            "ref_chrlist_path": ref_chrlist_path,
        }
        print(options)
        if data.marker_upload != "":  # 老师上传文件的string id
            upload_result = Dna("dna_gmap").find_one(collection="sg_upload_marker", query_dic={"_id": ObjectId(data.marker_upload)})
            marker_upload_path = upload_result['file_url']
            # options["marker_upload"] = os.path.join(base_path, marker_upload_path)
            options["marker_upload"] = Dna("dna_gmap").set_file_path(data.task_id, marker_upload_path, data.client)
        if data.generationid != "":   # 子代列表的id
            child_result = Dna("dna_gmap").find_one(collection="sg_child_list", query_dic={"_id": ObjectId(data.generationid)})
            options["child_list"] = ','.join(child_result['spcimen_ids'])
        self.set_sheet_data(name="dna_gmap.report.marker_filter", member_id=member_id, project_sn=project_sn,
                            task_id=data.task_id, main_table_name="Marker/" + main_table_name, options=options,
                            module_type="workflow", params=params, db_type="dna_gmap", analysis_name="Marker")
        task_info = super(MarkerFilterAction, self).POST()
        task_info['id'] = str(main_id)
        return json.dumps(task_info)


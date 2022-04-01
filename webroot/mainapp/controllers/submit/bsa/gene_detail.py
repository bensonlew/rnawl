# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# modified 2018.02.27

import re
import os
import web
import json
import datetime
from mainapp.libs.signature import check_sig
from mainapp.models.mongo.bsa import Bsa
from mainapp.controllers.project.bsa_controller import BsaController


class GeneDetailAction(BsaController):
    """
    BSA基因详情页接口, 暂时先改成 即时的分析
    lasted modified by hd 20180612
    """
    def __init__(self):
        super(GeneDetailAction, self).__init__(instant=True)

    @check_sig
    def POST(self):
        data = web.input()
        print data
        params = ["gene_id", "chrom", "start", "end", "type", "task_id", "task_type", "submit_location"]
        for param in params:
            if not hasattr(data, param):
                var = []
                var.append(param)
                info = {"success": False, "info": "缺少%s参数!" % param, "code": "C1200201", "variables": var }
                return json.dumps(info)
        if data.type not in ["detail", "seq"]:
            var =[]
            var.append(data.type)
            info = {"success": False, "info": "类型%s需为detail/seq!" % data.type, "code": "C1200202", "variables": var}
            return json.dumps(info)
        params_json = {
            "chrom": data.chrom,
            "gene_id": data.gene_id,
            "start": int(data.start),
            "end": int(data.end),
            "type": data.type,
            "task_type": int(data.task_type)
        }
        query_dic = {"task_id": data.task_id}
        result = Bsa().find_one_record(collection="sg_task", query_dic=query_dic)
        try:
            ref_fa_ = result["ref_path"]
            # index_file = result["index_path"]
            variant_file_ = result["variant_path"]
            project_sn = result["project_sn"]
            member_id = result["member_id"]
        except:
            info = {"success": False, "info": "sg_task表里没有ref_path、variant_path信息，请检查!",
                    "code": "C1200203", "variables": ""}
            return json.dumps(info)
        if data.type != "seq":
            params_json["submit_location"] = data.submit_location
        else:
            params_json["submit_location"] = data.submit_location
        ref_fa = Bsa().set_file_path(data.task_id, ref_fa_, data.client)
        variant_file = Bsa().set_file_path(data.task_id, variant_file_, data.client)
        params = json.dumps(params_json, sort_keys=True, separators=(',', ':'))
        if data.type == "seq":
            mongo_data = [
                ("name", "GeneSeq_" + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S"))),
                ("status", "start"),
                ("project_sn", project_sn),
                ("task_id", data.task_id),
                ("gene_id", data.gene_id),
                ("params", params),
                ("desc", "基因序列信息查找"),
                ("created_ts", datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
            ]
            main_id = Bsa().insert_main_table(collection="sg_gene_seq", data=mongo_data)
            Bsa().update_db_record(collection="sg_gene_seq", query_dict={"_id": main_id}, update_dict={"main_id": main_id})
            update_info = {str(main_id): "sg_gene_seq"}
        else:
            mongo_data = [
                ("name", "GeneDetail_" + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S"))),
                ("status", "start"),
                ("project_sn", project_sn),
                ("task_id", data.task_id),
                ("params", params),
                ("desc", "基因变异位点信息查找"),
                ("created_ts", datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
                ("gene_id", data.gene_id)
            ]
            main_id = Bsa().insert_main_table(collection="sg_gene_index", data=mongo_data)
            Bsa().update_db_record(collection="sg_gene_index", query_dict={"_id": main_id}, update_dict={"main_id": main_id})
            update_info = {str(main_id): "sg_gene_index"}
        options = {
            "gene_id": data.gene_id,
            "chrom": data.chrom,
            "start": int(data.start),
            "end": int(data.end),
            "type": data.type,
            "update_info": json.dumps(update_info),
            "main_id": str(main_id)
        }
        if data.type == "seq":
            options["ref_fa"] = ref_fa
        else:
            # options["index_file"] = index_file
            options["variant_file"] = variant_file
        sheet_data = self.set_sheet_data(name="bsa.gene_detail", member_id=member_id, project_sn=project_sn,
                                         task_id=data.task_id, main_table_name="gene_detail", options=options,
                                         module_type="tool", params=params)
        # print sheet_data
        task_info = super(GeneDetailAction, self).POST()
        task_info['id'] = str(main_id)
        return json.dumps(task_info)

# -*- coding: utf-8 -*-
# __author__ = 'HONGDONG'
# modified 20180422

import os
import re
import web
import json
import datetime
from mainapp.libs.signature import check_sig
from mainapp.models.mongo.dna import Dna
from mainapp.controllers.project.dna_controller import DnaController


class DrawCircosAction(DnaController):
    """
    画circos的接口
    """
    def __init__(self):
        super(DrawCircosAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        print data
        params = ["chrs", "task_id", "task_type", "submit_location"]
        for param in params:
            if not hasattr(data, param):
                var = []
                var.append(param)
                info = {"success": False, "info": "缺少%s参数!" % param, "code": "C3100501", "variables": var}
                return json.dumps(info)
        params_json = {
            "chrs": data.chrs,
            "task_type": int(data.task_type),
            "submit_location": data.submit_location
        }
        if not hasattr(data, "project_type"):
            db = "dna_wgs"
        elif data.project_type == "dna_gmap":
            db = "dna_gmap"
        else:
            db = data.project_type
        result = Dna(db).find_one(collection="sg_task", query_dic={"task_id": data.task_id})
        # noinspection PyBroadException
        try:
            cnv_anno_path = result['cnv_anno_path']
            sv_anno_path = result['sv_anno_path']
            snp_anno = result['snp_anno_vcf']  # snp.anno.primary.vcf
            indel_anno = result['indel_anno_vcf']  # indel.anno.primary.vcf
            project_sn = result["project_sn"]
            member_id = result["member_id"]
            genome_version_id = result['genome_version_id']
        except:
            info = {"success": False, "info": "sg_task表里没有project_sn, member_id, species_version_id信息，请检查!",
                    "code": "C3100502", "variables": ""}
            return json.dumps(info)
        if cnv_anno_path.startswith("rerewrweset"):
            base_path = '/mnt/ilustre/data/' if str(data.client) == 'client01' else "/mnt/ilustre/tsanger-data/"
            cnv_anno_path_ = os.path.join(base_path, cnv_anno_path)
            sv_anno_path_ = os.path.join(base_path, sv_anno_path)
            snp_anno_ = os.path.join(base_path, snp_anno)
            indel_anno_ = os.path.join(base_path, indel_anno)
        elif cnv_anno_path.startswith("//"):
            region = result['region']
            cnv_anno_path_ = ":".join([region.rstrip(':'), cnv_anno_path])
            sv_anno_path_ = ":".join([region.rstrip(':'), sv_anno_path])
            snp_anno_ = ":".join([region.rstrip(':'), snp_anno])
            indel_anno_ = ":".join([region.rstrip(':'), indel_anno])
        elif cnv_anno_path.startswith("/mnt") or re.match(".*://.*", cnv_anno_path):
            cnv_anno_path_ = cnv_anno_path
            sv_anno_path_ = sv_anno_path
            snp_anno_ = snp_anno
            indel_anno_ = indel_anno
        else:
            info = {"success": False, "info": "存入的文件%s路径格式不正确！", "variables":cnv_anno_path, "code" : "C3100504"}
            return json.dumps(info)
        params = json.dumps(params_json, sort_keys=True, separators=(',', ':'))
        main_table_name = 'Circos_' + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        mongo_data = [
            ("name", main_table_name),
            ("status", "start"),
            ("project_sn", project_sn),
            ("task_id", data.task_id),
            ("params", params),
            ("desc", "Circos图主表！"),
            ("created_ts", datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        ]
        main_id = Dna(db).insert_main_table(collection="sg_circos", data=mongo_data)
        Dna(db).update_db_record(collection="sg_circos", query_dict={"_id": main_id},
                                        update_dict={"main_id": main_id})
        update_info = {str(main_id): "sg_circos"}
        options = {
            "windows": 100000,
            "snp": snp_anno_,
            "indel": indel_anno_,
            "sv_path": sv_anno_path_,
            "cnv_path": cnv_anno_path_,
            "chrs": data.chrs,
            "update_info": json.dumps(update_info),
            "main_id": str(main_id),
            "genome_version_id": str(genome_version_id)
        }
        if hasattr(data, "project_type"):
            options["project_type"] = data.project_type
        self.set_sheet_data(name="wgs.report.draw_circos", member_id=member_id, project_sn=project_sn,
                            task_id=data.task_id, main_table_name="circos/" + main_table_name,
                            options=options, params=params, db_type=db)
        task_info = super(DrawCircosAction, self).POST()
        task_info['id'] = str(main_id)
        return json.dumps(task_info)

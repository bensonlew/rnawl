# -*- coding: utf-8 -*-
# __author__ = 'HONGDONG'
# modified 20180425

import re
import os
import web
import json
import datetime
from biocluster.config import Config
from mainapp.libs.signature import check_sig
from mainapp.models.mongo.dna import Dna
from mainapp.controllers.project.dna_controller import DnaController


class AssemblyAction(DnaController):
    """
    画circos的接口
    """
    def __init__(self):
        super(AssemblyAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        print data
        params = ["group_dict", "group_id", "poss", "unmapping", "task_id", "task_type", "submit_location"]
        for param in params:
            if not hasattr(data, param):
                var = []
                var.append(param)
                info = {"success": False, "info": "缺少%s参数!" % param, "code": "C3100101", "variables": var}
                return json.dumps(info)
        if data.unmapping not in ["true", "false"]:
            var = []
            var.append(data.unmapping)
            info = {"success": False, "info": "参数 %s 不合法, 必须为true or false！" % data.unmapping,
                    "code": "C3100102", "variables": var}
            return json.dumps(info)
        params_json = {
            "group_dict": json.loads(data.group_dict),
            "group_id": data.group_id,
            "poss": data.poss,
            "unmapping": data.unmapping,
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
            bam_path = result['bam_path']
            project_sn = result["project_sn"]
            member_id = result["member_id"]
            genome_version_id = result['genome_version_id']
        except:
            info = {"success": False, "info": "sg_task表里没有project_sn, member_id信息，请检查!",
                    "code": "C3100103", "variables": ""}
            return json.dumps(info)
        is_bucket = "false"
        if bam_path.startswith("rerewrweset"):
            is_bucket = "true"
            base_path = 'http://bcl.sanger.com/data/' if str(data.client) == 'client01' \
                else "http://bcl.tsanger.com/data/"
            bam_path_ = base_path + bam_path.lstrip("/")
        elif bam_path.startswith("//"):
            is_bucket = "true"
            region = result['region']
            bam_path_ = ":".join([region.rstrip(':'), bam_path])
        elif bam_path.startswith("/mnt"):
            bam_path_ = bam_path
        elif re.match(".*://.*", bam_path):
            is_bucket = "true"
            bam_path_ = bam_path
        else:
            info = {"success": False, "info": "存入的文件%s路径格式不正确！", "variables":bam_path, "code" : "C3100104"}
            return json.dumps(info)
        params = json.dumps(params_json, sort_keys=True, separators=(',', ':'))
        main_table_name = 'Assembly_statistic' + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        mongo_data = [
            ("name", main_table_name),
            ("status", "start"),
            ("project_sn", project_sn),
            ("task_id", data.task_id),
            ("params", params),
            ("desc", "Assembly主表！"),
            ("created_ts", datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        ]
        main_id = Dna(db).insert_main_table(collection="sg_assembly", data=mongo_data)
        Dna(db).update_db_record(collection="sg_assembly", query_dict={"_id": main_id},
                                 update_dict={"main_id": main_id})
        update_info = {str(main_id): "sg_assembly"}
        options = {
            "samples": ','.join(json.loads(data.group_dict)['all']),
            "poss": data.poss,
            "unmapping": data.unmapping,
            "bam_path": bam_path_,
            "update_info": json.dumps(update_info),
            "main_id": str(main_id),
            "species_version_id": str(genome_version_id),
            "is_bucket": is_bucket
        }
        if hasattr(data, "project_type"):
            options["project_type"] = data.project_type
        self.set_sheet_data(name="wgs.report.assembly", member_id=member_id, project_sn=project_sn,
                            task_id=data.task_id, main_table_name="assembly/" + main_table_name,
                            options=options, params=params, db_type=db)
        task_info = super(AssemblyAction, self).POST()
        task_info['id'] = str(main_id)
        return json.dumps(task_info)

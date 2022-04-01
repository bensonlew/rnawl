# -*- coding: utf-8 -*-
# __author__ = 'HONGDONG'
# modified 20190308

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
    有参变异检测V2版中的局部组装的接口
    """
    def __init__(self):
        super(AssemblyAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        print data
        params = ["samples", "poss", "posname", "unmapping", "task_id", "task_type", "submit_location",
                  "chongmingming_result"]
        for param in params:
            if not hasattr(data, param):
                info = {"success": False, "info": "缺少%s参数!" % param}
                return json.dumps(info)
        if not data.posname:
            info = {"success": False, "info": "请对每个区域进行命名！"}
            return json.dumps(info)
        if data.unmapping not in ["true", "false"]:
            info = {"success": False, "info": "参数 %s 不合法, 必须为true or false！" % data.unmapping}
            return json.dumps(info)
        if len(data.poss.split('|')) != len(data.posname.split('|')):
            info = {"success": False, "info": "区域的个数与区域名字的个数不对应！"}
            return json.dumps(info)
        params_json = {
            # "group_dict": json.loads(data.group_dict),
            # "group_id": data.group_id,
            "samples": data.samples,
            "posname": data.posname,
            "poss": data.poss,
            "unmapping": data.unmapping,
            "task_type": int(data.task_type),
            "submit_location": data.submit_location,
            "task_id": data.task_id,
            "chongmingming_result": data.chongmingming_result
        }
        if not hasattr(data, "project_type"):
            db = "dna_wgs_v2"
        elif data.project_type == "dna_gmap":
            db = "dna_gmap"
        elif data.project_type == "dna_wgs":
            db = "dna_wgs"
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
            info = {"success": False, "info": "sg_task表里没有project_sn, member_id信息，请检查!"}
            return json.dumps(info)
        bam_path_ = Dna(db).set_file_path(data.task_id, bam_path, data.client)
        is_bucket = False if not bam_path_.startswith('s3') else True
        params = json.dumps(params_json, sort_keys=True, separators=(',', ':'))
        # main_table_name = 'Assembly_statistic_' + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        main_table_name = Dna(db).set_main_table_name("Assembly_statistic", data.chongmingming_result)
        mongo_data = [
            ("name", main_table_name),
            ("status", "start"),
            ("project_sn", project_sn),
            ("task_id", data.task_id),
            ("params", params),
            ("desc", "Assembly主表！"),
            ('member_id', member_id),
            ("created_ts", datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        ]
        main_id = Dna(db).insert_main_table(collection="sg_assembly", data=mongo_data)
        Dna(db).update_db_record(collection="sg_assembly", query_dict={"_id": main_id},
                                 update_dict={"main_id": main_id})
        update_info = {str(main_id): "sg_assembly"}
        options = {
            # "samples": ','.join(json.loads(data.group_dict)['all']),
            "samples": data.samples,
            "poss": data.poss,
            "unmapping": data.unmapping,
            "bam_path": bam_path_,
            "update_info": json.dumps(update_info),
            "main_id": str(main_id),
            "species_version_id": str(genome_version_id),
            "is_bucket": is_bucket,
            "posname": data.posname
        }
        if hasattr(data, "project_type"):
            options["project_type"] = data.project_type
        self.set_sheet_data(name="wgs_v2.report.assembly", member_id=member_id, project_sn=project_sn,
                            task_id=data.task_id, main_table_name="assembly/" + main_table_name,
                            options=options, params=params, db_type=db)
        task_info = super(AssemblyAction, self).POST()
        task_info['id'] = str(main_id)
        return json.dumps(task_info)

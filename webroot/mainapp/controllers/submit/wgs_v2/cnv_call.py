# -*- coding: utf-8 -*-
# __author__ = 'Liuwentian'
# modified 20190308

import os
import web
import json
import datetime
from mainapp.libs.signature import check_sig
from mainapp.models.mongo.dna import Dna
from mainapp.controllers.project.dna_controller import DnaController
from bson.objectid import ObjectId


class CnvCallAction(DnaController):
    """
     cnv统计接口
    """
    def __init__(self):
        super(CnvCallAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        print data
        params = ["task_id", "task_type", "submit_location", "chongmingming_result"]
        for param in params:
            if not hasattr(data, param):
                info = {"success": False, "info": "缺少%s参数!" % param}
                return json.dumps(info)
        db = data.project_type = "dna_wgs_v2"
        result = Dna(db).find_one(collection="sg_task", query_dic={"task_id": data.task_id})
        try:
            project_sn = result["project_sn"]
            genome_version_id = result["genome_version_id"]
            bam_path = result["bam_path"]
            member_id = result["member_id"]
        except:
            info = {"success": False, "info": "sg_task表里没有project_sn, genome_version_id,bam_path信息，请检查!"}
            return json.dumps(info)
        ref_result = Dna(db).find_one(collection="sg_species_version", query_dic={"_id": genome_version_id})
        try:
            ref_gff = ref_result['gff']
        except:
            info = {"success": False, "info": "sg_species_version表里没有ref、snpeff_path信息，请检查!"}
            return json.dumps(info)
        params_json = {
            "task_type": int(data.task_type),
            "submit_location": data.submit_location,
            "task_id": data.task_id,
            "chongmingming_result": data.chongmingming_result
        }
        params = json.dumps(params_json, sort_keys=True, separators=(',', ':'))
        main_table_name = Dna(db).set_main_table_name("Cnv_call", data.chongmingming_result)
        mongo_data = [
            ("name", main_table_name),
            ("status", "start"),
            ("project_sn", project_sn),
            ("task_id", data.task_id),
            ("params", params),
            ("desc", "Cnv_call主表！"),
            ("created_ts", datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        ]
        main_id = Dna(db).insert_main_table(collection="sg_cnv_call", data=mongo_data)
        Dna(db).update_db_record(collection="sg_cnv_call", query_dict={"_id": main_id},
                                 update_dict={"main_id": main_id})
        update_info = {str(main_id): "sg_cnv_call"}
        options = {
            "ref_gff": ref_gff,
            "bam_list": bam_path,
            "update_info": json.dumps(update_info),
            "main_id": str(main_id),
            "task_id": data.task_id
        }
        main_table_name = "Cnv_call_" + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3])
        self.set_sheet_data(name="wgs_v2.report.cnv_call", member_id=member_id, project_sn=project_sn,
                            task_id=data.task_id,
                            main_table_name=main_table_name, options=options, module_type="workflow", params=params,
                            target_output=True, db_type=db)
        task_info = super(CnvCallAction, self).POST()
        task_info['id'] = str(main_id)
        return json.dumps(task_info)

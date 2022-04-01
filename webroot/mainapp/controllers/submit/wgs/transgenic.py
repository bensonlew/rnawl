# -*- coding: utf-8 -*-
# __author__ = 'Liuwentian'
# modified 20180518

import os
import web
import json
import datetime
from mainapp.libs.signature import check_sig
from mainapp.models.mongo.dna import Dna
from mainapp.controllers.project.dna_controller import DnaController
from bson.objectid import ObjectId

class TransgenicAction(DnaController):
    """
     转基因接口
    """
    def __init__(self):
        super(TransgenicAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        print data
        params = ["group_dict", "group_id", "transgenic_file_id", "task_id", "task_type", "submit_location"]
        for param in params:
            if not hasattr(data, param):
                var = []
                var.append(param)
                info = {"success": False, "info": "缺少%s参数!" % param, "code": "C3101401", "variables": var}
                return json.dumps(info)
        params_json = {
            "group_dict": json.loads(data.group_dict),
            "group_id": data.group_id,
            "transgenic_file_id": data.transgenic_file_id,
            "task_type": int(data.task_type),
            "submit_location": data.submit_location
        }
        result = Dna("dna_wgs").find_one(collection="sg_task", query_dic={"task_id": data.task_id})
        try:
            project_sn = result["project_sn"]
            genome_version_id = result['genome_version_id']
            clean_fastq_path = result["clean_fastq_path"]
            member_id = result["member_id"]
        except:
            info = {"success": False, "info": "sg_task表里没有genome_version_id, clean_fastq_path信息，请检查!",
                    "code": "C3101402", "variables": ""}
            return json.dumps(info)
        try:
            record_dict = Dna("dna_wgs").find_one(collection="sg_upload_transgenic", query_dic={"_id": ObjectId(data.transgenic_file_id)})
            file_url = record_dict["file_url"]
            if str(data.client) == 'client01':
                base_path = '/mnt/ilustre/data/'
            else:
                base_path = "/mnt/ilustre/tsanger-data/"
            insert_path = os.path.join(base_path, file_url)
        except:
            info = {"success": False, "info": "sg_upload_transgenic表里没有transgenic_file_id信息，请检查!",
                    "code": "C3101403", "variables": ""}
            return json.dumps(info)
        params = json.dumps(params_json, sort_keys=True, separators=(',', ':'))
        main_table_name = 'Transgene_statistic' + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        mongo_data = [
            ("name", main_table_name),
            ("status", "start"),
            ("project_sn", project_sn),
            ("task_id", data.task_id),
            ("params", params),
            ("desc", "Transgene主表！"),
            ("created_ts", datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        ]
        main_id = Dna("dna_wgs").insert_main_table(collection="sg_transgenic", data=mongo_data)
        Dna("dna_wgs").update_db_record(collection="sg_transgenic", query_dict={"_id": main_id},
                                        update_dict={"main_id": main_id})
        update_info = {str(main_id): "sg_transgenic"}
        options = {
            "samples": json.dumps(json.loads(data.group_dict)['all']),
            "clean_path": clean_fastq_path,
            "insert_seq": insert_path,
            "genome_version_id": str(genome_version_id),
            "update_info": json.dumps(update_info),
            "main_id": str(main_id)
        }
        self.set_sheet_data(name="wgs.report.transgenic", member_id=member_id, project_sn=project_sn,
                            task_id=data.task_id, main_table_name="transgenic/" + main_table_name,
                            options=options, params=params)
        task_info = super(TransgenicAction, self).POST()
        task_info['id'] = str(main_id)
        return json.dumps(task_info)
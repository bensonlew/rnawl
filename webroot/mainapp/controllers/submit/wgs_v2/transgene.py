# -*- coding: utf-8 -*-
# __author__ = 'Liuwentian'
# modified 20190312

import os
import web
import json
import datetime
from mainapp.libs.signature import check_sig
from mainapp.models.mongo.dna import Dna
from mainapp.controllers.project.dna_controller import DnaController
from bson.objectid import ObjectId


class TransgeneAction(DnaController):
    """
     转基因
    """
    def __init__(self):
        super(TransgeneAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        print data
        params = ["sample", "insert_id", "task_id", "task_type", "submit_location", "chongmingming_result"]
        for param in params:
            if not hasattr(data, param):
                info = {"success": False, "info": "缺少%s参数!" % param}
                return json.dumps(info)
        db = data.project_type = "dna_wgs_v2"
        result = Dna(db).find_one(collection="sg_task", query_dic={"task_id": data.task_id})
        try:
            project_sn = result["project_sn"]
            genome_version_id = result["genome_version_id"]
            # clean_fastq_path = result["clean_fastq_path"]
            raw_data_path = result["raw_data_path"]
            member_id = result["member_id"]
        except:
            info = {"success": False, "info": "sg_task表里没有project_sn, genome_version_id,raw_data_path,member_id信息，请检查!"}
            return json.dumps(info)
        rawdata_result = Dna(db).find_one(collection="sg_specimen_other", query_dic={"task_id": data.task_id,
                                                                                     "analysis_name": data.sample})
        try:
            file_name = rawdata_result["file_name"]
        except:
            info = {"success": False, "info": "sg_specimen_other表里没有file_name信息，请检查!"}
            return json.dumps(info)
        file_name_list = file_name.strip().split(",")
        if len(file_name_list) == 1:
            rawdata_path_1 = os.path.join(raw_data_path, file_name_list[0])
            rawdata_path_2 = ""
        else:
            rawdata_path_1 = os.path.join(raw_data_path, file_name_list[0])
            rawdata_path_2 = os.path.join(raw_data_path, file_name_list[1])
        ref_result = Dna(db).find_one(collection="sg_species_version", query_dic={"_id": genome_version_id})
        try:
            ref_fa = ref_result['ref']
        except:
            info = {"success": False, "info": "sg_species_version表里没有ref信息，请检查!"}
            return json.dumps(info)
        insert_result = Dna(db).find_one(collection="sg_feature_file", query_dic={"_id": ObjectId(data.insert_id)})
        try:
            insertion_id = insert_result['insertion_id']
            sequence = insert_result['sequence']
        except:
            info = {"success": False, "info": "sg_feature_file表里没有insertion_id,sequence信息，请检查!"}
            return json.dumps(info)
        main_table_name = Dna(db).set_main_table_name("Transgene", data.chongmingming_result)
        params_json = {
            "sample": data.sample,
            "insert_id": data.insert_id,
            "task_type": int(data.task_type),
            "submit_location": data.submit_location,
            "task_id": data.task_id,
            "chongmingming_result": data.chongmingming_result
        }
        params = json.dumps(params_json, sort_keys=True, separators=(',', ':'))
        mongo_data = [
            ("name", main_table_name),
            ("status", "start"),
            ("project_sn", project_sn),
            ("task_id", data.task_id),
            ("params", params),
            ("desc", "Transgene主表！"),
            ("created_ts", datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ("sample", data.sample)
        ]
        main_id = Dna(db).insert_main_table(collection="sg_transgenosis", data=mongo_data)
        Dna(db).update_db_record(collection="sg_transgenosis", query_dict={"_id": main_id},
                                 update_dict={"main_id": main_id})
        update_info = {str(main_id): "sg_transgenosis"}
        options = {
            "ref_fa": ref_fa,
            "insertion_id": insertion_id,
            "sequence": sequence,
            "rawdata_path_1": rawdata_path_1,
            "sample": data.sample,
            "update_info": json.dumps(update_info),
            "main_id": str(main_id),
            "task_id": data.task_id
        }
        if rawdata_path_2 != "":
            options["rawdata_path_2"] = rawdata_path_2
        main_table_name = "Transgene_" + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3])
        self.set_sheet_data(name="wgs_v2.report.transgene", member_id=member_id, project_sn=project_sn,
                            task_id=data.task_id,
                            main_table_name=main_table_name, options=options, module_type="workflow", params=params,
                            db_type=db)
        task_info = super(TransgeneAction, self).POST()
        task_info['id'] = str(main_id)
        return json.dumps(task_info)

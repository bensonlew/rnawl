# -*- coding: utf-8 -*-
# __author__ = 'HONGDONG'
# modified 20180515

import os
import web
import json
import datetime
from bson.objectid import ObjectId
from mainapp.libs.signature import check_sig
from mainapp.models.mongo.dna import Dna
from mainapp.controllers.project.dna_controller import DnaController


class BlastSeqAction(DnaController):
    """
    局部组装与转基因用到的blast与序列获取接口
    """
    def __init__(self):
        super(BlastSeqAction, self).__init__(instant=True)

    @check_sig
    def POST(self):
        data = web.input()
        print data
        params = ["main_id", "sample_id", "scaffold_id", "analysis_type", "blast_type", "task_id",
                  "task_type", "submit_location"]
        for param in params:
            if not hasattr(data, param):
                var = []
                var.append(param)
                info = {"success": False, "info": "缺少%s参数!" % param, "code": "C3100201", "variables": var}
                return json.dumps(info)
        params_json = {
            "sample_id": data.sample_id,
            "scaffold_id": data.scaffold_id,
            "analysis_type": data.analysis_type,
            "blast_type": data.blast_type,
            "task_type": int(data.task_type),
            "main_id": data.main_id,
            "submit_location": data.submit_location
        }
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
            genome_version_id = result['genome_version_id']
        except:
            info = {"success": False, "info": "sg_task表里没有project_sn, member_id, species_version_id信息，请检查!",
                    "code": "C3100202", "variables": ""}
            return json.dumps(info)
        if data.analysis_type == 'blast':
            if data.blast_type == "ref":
                genome_info = Dna(db).find_one(collection="sg_species_version",
                                                      query_dic={"_id": genome_version_id})
                try:
                    ref_db_path = genome_info['ssr_path']
                except:
                    var = []
                    var.append(genome_version_id)
                    info = {"success": False,
                            "info": "sg_species_version表里没有%s信息，请检查!" % (genome_version_id),
                            "code": "C3100203", "variables": var}
                    return json.dumps(info)
                if data.client == "client03":
                    ref_db = os.path.join("/mnt/ilustre/users/sanger-dev/app/database/dna_wgs_geneome/",
                                          "{}/makedbblast/ref".format(ref_db_path.strip('/')))
                else:
                    ref_db = os.path.join("/mnt/ilustre/users/sanger/app/database/dna_wgs_geneome/",
                                          "{}/makedbblast/ref".format(ref_db_path.strip('/')))
            else:
                tran_info = Dna(db).find_one(collection="sg_transgene",
                                                    query_dic={"_id": ObjectId(data.main_id)})
                try:
                    ref_db = tran_info['insert_path']
                except:
                    info = {"success": False, "info": "sg_transgene表里没有insert_path字段信息，请检查!",
                            "code": "C3100204", "variables": ""}
                    return json.dumps(info)
        else:
            ref_db = 'no use'
        if data.submit_location == "blastseq":  # 这个要重新定义与前端一致
            collection = "sg_assembly"
        else:
            collection = "sg_transgene"
        seq_info = Dna(db).find_one(collection=collection, query_dic={"_id": ObjectId(data.main_id)})
        # noinspection PyBroadException
        try:
            insert_seq = ''
            # seq_path = os.path.join('s3://', seq_info['seq_path'])
            seq_path = seq_info['seq_path']
            if collection == "sg_transgene":
                insert_seq = seq_info['insert_seq']
        except:
            var = []
            var.append(collection)
            var.append(data.main_id)
            info = {"success": False, "info": "%s表里没有%s信息，请检查!" % (collection, data.main_id),
                    "code": "C3100205", "variables": var}
            return json.dumps(info)
        seq_path = Dna(db).set_file_path(data.task_id, seq_path, data.client)
        params = json.dumps(params_json, sort_keys=True, separators=(',', ':'))
        main_table_name = "{}_".format(data.analysis_type.capitalize()) + \
                          datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        mongo_data = [
            ("name", main_table_name),
            ("status", "start"),
            ("project_sn", project_sn),
            ("task_id", data.task_id),
            ("params", params),
            ("types", data.analysis_type),
            ("desc", "--"),
            ("created_ts", datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        ]
        main_id = Dna(db).insert_main_table(collection="sg_blast_seq", data=mongo_data)
        update_info = {str(main_id): "sg_blast_seq"}
        options = {
            "sample_id": data.sample_id,
            "scaffold_id": data.scaffold_id,
            "seq_path": seq_path.rstrip('/') + '/',
            "types": data.analysis_type,
            "ref_db": ref_db,
            "blast_type": data.blast_type,
            "update_info": json.dumps(update_info),
            "main_id": str(main_id)
        }
        if hasattr(data, "project_type"):
            options["project_type"] = data.project_type
        if data.analysis_type == 'blast' and data.blast_type == 'insert':
            options.update({"insert_seq": insert_seq})
        self.set_sheet_data(name="wgs.report.blast_seq", member_id=member_id, project_sn=project_sn,
                            task_id=data.task_id, main_table_name="blast_seq/" + main_table_name,
                            options=options, params=params, db_type=db)
        task_info = super(BlastSeqAction, self).POST()
        task_info['id'] = str(main_id)
        return json.dumps(task_info)

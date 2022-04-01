# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# modified 2018.05.07

import os
import re
import web
import json
import datetime
from mainapp.libs.signature import check_sig
from mainapp.models.mongo.dna import Dna
from mainapp.controllers.project.dna_controller import DnaController


class GeneStructureAction(DnaController):
    """
    基因详情页-基因结构接口
    """
    def __init__(self):
        super(GeneStructureAction, self).__init__(instant=False)

    def POST(self):
        data = web.input()
        print data
        params = ["gene_id", "chrom", "start", "end", "task_id", "task_type", "submit_location"]
        for param in params:
            if not hasattr(data, param):
                var = []
                var.append(param)
                info = {"success": False, "info": "缺少%s参数!" % param, "code": "C3100701", "variables": var}
                return json.dumps(info)
        if not data.gene_id:
            info = {"success": False, "info": "参数gene_id不能为空!", "code": "C3100702", "variables": ""}
            return json.dumps(info)
        if not data.chrom:
            info = {"success": False, "info": "参数gene_id不能为空!", "code": "C3100703", "variables": ""}
            return json.dumps(info)
        if int(data.start) >= int(data.end):
            info = {"success": False, "info": "参数start不能小于等于end!", "code": "C3100704", "variables": ""}
            return json.dumps(info)
        if not hasattr(data, "project_type"):
            db = "dna_wgs"
        elif data.project_type == "dna_gmap":
            db = "dna_gmap"
        else:
            db = data.project_type
        result = Dna(db).find_one(collection="sg_task", query_dic={"task_id": data.task_id})
        # noinspection PyBroadException
        try:
            pop_final_vcf = result['pop_final_vcf']
            project_sn = result["project_sn"]
            member_id = result["member_id"]
            genome_version_id = result["genome_version_id"]
        except:
            info = {"success": False, "info": "sg_task表里没有project_sn, member_id,pop_final_vcf信息，请检查!",
                    "code": "C3100705", "variables": ""}
            return json.dumps(info)
        pop_final_vcf = Dna(db).set_file_path(data.task_id, pop_final_vcf, data.client)
        ref_result = Dna(db).find_one(collection="sg_species_version", query_dic={"_id": genome_version_id})
        # noinspection PyBroadException
        try:
            ref_fa = ref_result['ref']
            ref_gff = ref_result["gff"]
        except:
            info = {"success": False, "info": "sg_species_version表里没有ref、gff信息，请检查!",
                    "code": "C3100707", "variables": ""}
            return json.dumps(info)
        params_json = {
            "gene_id": data.gene_id,
            "chrom": data.chrom,
            "start": int(data.start),
            "end": int(data.end),
            "task_type": int(data.task_type),
            "submit_location": data.submit_location
        }
        params = json.dumps(params_json, sort_keys=True, separators=(',', ':'))
        mongo_data = [
            ("name", "GeneStructure_" + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S"))),
            ("status", "start"),
            ("project_sn", project_sn),
            ("task_id", data.task_id),
            ("params", params),
            ("gene_id", data.gene_id),
            ("desc", "基因详情页-基因结构"),
            ("created_ts", datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        ]
        main_id = Dna(db).insert_main_table(collection="sg_gene_structure", data=mongo_data)
        Dna(db).update_db_record(collection="sg_gene_structure", query_dict={"_id": main_id},
                                 update_dict={"main_id": main_id})
        update_info = {str(main_id): "sg_gene_structure"}
        options = {
            "ref_fa": ref_fa,
            "ref_gff": ref_gff,
            "pop_final_vcf": pop_final_vcf,
            "chrom": data.chrom,
            "location": data.start + "-" + data.end,
            "gene_id": data.gene_id,
            "update_info": json.dumps(update_info),
            "main_id": str(main_id)
        }
        if hasattr(data, "project_type"):
            options["project_type"] = data.project_type
        main_table_name = "GeneDetail/GeneStructure_" + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3])
        self.set_sheet_data(name="wgs.gene_structure", member_id=member_id, project_sn=project_sn, task_id=data.task_id,
                            main_table_name=main_table_name, options=options, module_type="tool", params=params,
                            db_type=db)
        task_info = super(GeneStructureAction, self).POST()
        task_info['id'] = str(main_id)
        return json.dumps(task_info)

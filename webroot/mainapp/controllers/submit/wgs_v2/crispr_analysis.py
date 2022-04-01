# -*- coding: utf-8 -*-
# __author__ = 'Binbin Zhao
# modified 20190403


import web
import json
import datetime
import os
from bson.objectid import ObjectId
from types import StringTypes
from mainapp.libs.signature import check_sig
from mainapp.models.mongo.dna import Dna
from mainapp.controllers.project.dna_controller import DnaController


class CrisprAnalysisAction(DnaController):
    """
    wgs_v2 染色体分布覆盖分布图接口
    """

    def __init__(self):
        super(CrisprAnalysisAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        print data
        params = ["task_id", "task_type", "wt", "bam_list", "region", "submit_location", "sg_rna",
                  "chongmingming_result", "gene_name"]
        for param in params:
            if not hasattr(data, param):
                info = {"success": False, "info": "缺少%s参数!" % param}
                return json.dumps(info)
        task_result = Dna("dna_wgs_v2").find_one(collection="sg_task", query_dic={"task_id": data.task_id})
        if not task_result:
            info = {"success": False,
                    "info": "sg_task表里没有task_id: %s对应的信息，请检查!" % data.task_id}
            return json.dumps(info)
        project_sn = task_result["project_sn"]
        member_id = task_result["member_id"]
        bam_path = task_result["bam_path"]
        fa_path = Dna("dna_wgs_v2").find_one(collection="sg_species_version",
                                             query_dic={"_id": task_result["genome_version_id"]})
        if not fa_path:
            info = {"success": False,
                    "info": "sg_species_version表里没有_id: %s对应的信息，请检查!" % data.file_id}
            return json.dumps(info)
        else:
            try:
                file_path = fa_path['ref']
            except:
                info = {"success": False,
                        "info": "dep_path中没有dep_path字段，请核查！"}
                return json.dumps(info)
        # file_path = Dna("dna_wgs_v2").set_file_path(data.task_id, file_path, data.client)
        params_json = {
            "task_type": int(data.task_type),
            "wt": data.wt,
            "bam_list": json.loads(data.bam_list),
            "region": data.region,
            "sg_rna": data.sg_rna,
            "submit_location": data.submit_location,
            "chongmingming_result": data.chongmingming_result,
            "gene_name": data.gene_name,
            "task_id": data.task_id
        }
        params = json.dumps(params_json, sort_keys=True, separators=(',', ':'))     # 返回json
        main_table_name = 'Crispr_Analysis_' + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")
        mongo_data = [
            ("name", main_table_name),
            ("project_sn", project_sn),
            ("task_id", data.task_id),
            ("created_ts", datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ("status", "start"),
            ("params", params),
            ("desc", "基因编辑CRISPR主表"),
        ]
        main_id = Dna("dna_wgs_v2").insert_main_table(collection="sg_crispr_analysis", data=mongo_data)
        Dna("dna_wgs_v2").update_db_record(collection="sg_crispr_analysis", query_dict={"_id": main_id},
                                           update_dict={"main_id": main_id})
        update_info = {str(main_id): "sg_crispr_analysis"}
        options = {
            "project_sn": project_sn,
            "main_id": str(main_id),
            "task_id": data.task_id,
            "update_info": json.dumps(update_info),     # 给work_flow用
            "wt": data.wt,
            "bam_list": data.bam_list,
            "region": data.region,
            "sg_rna": data.sg_rna,
            "gene_name": data.gene_name,
            "ref_fa": file_path,
            "bam_path": bam_path
        }
        self.set_sheet_data(name="wgs_v2.crispr_analysis", member_id=member_id, project_sn=project_sn,
                            task_id=data.task_id, main_table_name="CrisprAnalysis/" + main_table_name, options=options,
                            module_type="module", params=params, db_type="dna_wgs_v2", analysis_name="CrisprAnalysis")
        task_info = super(CrisprAnalysisAction, self).POST()
        task_info['id'] = str(main_id)
        return json.dumps(task_info)

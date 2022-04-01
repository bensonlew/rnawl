# -*- coding: utf-8 -*-
# __author__ = 'wentianliu'
# modified 20180919
import web
import json
import datetime
import os
from bson.objectid import ObjectId
from types import StringTypes
from mainapp.libs.signature import check_sig
from mainapp.models.mongo.dna import Dna
from mainapp.controllers.project.dna_controller import DnaController
from biocluster.file import getsize, exists, list_dir
from biocluster.file import get_reader, get_top_lines, get_last_lines


class GeneSampleSeqAction(DnaController):
    """
    基因详情页-样本序列接口
    """
    def __init__(self):
        super(GeneSampleSeqAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        print data
        params = ["gene_id", "group_dict", "task_id", "task_type", "submit_location"]
        for param in params:
            if not hasattr(data, param):
                info = {"success": False, "info": "缺少%s参数!" % param}
                return json.dumps(info)
        group_dict_ = json.loads(data.group_dict)
        specimen_ids = []
        for k in group_dict_.keys():
            for s in group_dict_[k]:
                specimen_ids.append(s)
        specimen_ids = list(set(specimen_ids))
        # data.project_type = "dna_evolution"
        if not hasattr(data, "project_type"):
            db = "dna_wgs"
        elif data.project_type == "dna_gmap":
            db = "dna_gmap"
        elif data.project_type == "dna_wgs_v2":
            db = "dna_wgs_v2"
        else:
            db = data.project_type
        result = Dna(db).find_one(collection="sg_task", query_dic={"task_id": data.task_id})
        print db
        try:
            pop_final_vcf = result['pop_final_vcf']
            project_sn = result["project_sn"]
            member_id = result["member_id"]
            genome_version_id = result["genome_version_id"]
            if db == "dna_wgs_v2":
                is_wgs_result = "dna_wgs_v2"
            else:
                is_wgs_result = result["is_wgs_result"]
        except:
            info = {"success": False, "info": "sg_task表里没有project_sn, member_id,pop_final_vcf信息，请检查!"}
            return json.dumps(info)
        if str(data.client) == 'client01':
            base_path = '/mnt/ilustre/data/'
        else:
            base_path = "/mnt/ilustre/tsanger-data/"
        if pop_final_vcf.startswith("rerewrweset"):
            pop_final_vcf = os.path.join(base_path, pop_final_vcf)
        if not exists(pop_final_vcf):
        # if not os.path.exists(pop_final_vcf):
            info = {"success": False, "info": "sg_task表里的pop_final_vcf:{}文件不存在，请检查!".format(pop_final_vcf)}
            return json.dumps(info)
        if is_wgs_result == "yes":
            mongo_dna = "dna_wgs"
        elif is_wgs_result == "dna_wgs_v2":
            mongo_dna = "dna_wgs_v2"
        else:
            mongo_dna = "dna_evolution"
        ref_result = Dna(mongo_dna).find_one(collection="sg_species_version", query_dic={"_id": genome_version_id})
        try:
            ref_fa = ref_result['ref']
            ref_gff = ref_result["gff"]
            ssr_path = ref_result["ssr_path"]
        except:
            info = {"success": False, "info": "sg_species_version表里没有ref、gff信息，请检查!"}
            return json.dumps(info)
        params_json = {
            "gene_id": data.gene_id,
            "group_id": data.group_id,
            "group_dict": group_dict_,
            "task_type": int(data.task_type),
            "submit_location": data.submit_location,
            "task_id": data.task_id,
            "anno_id": data.anno_id,
            "gene_name": data.gene_name
        }
        params = json.dumps(params_json, sort_keys=True, separators=(',', ':'))
        mongo_data = [
            ("name", "SamplesSeq_" + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S"))),
            ("status", "start"),
            ("project_sn", project_sn),
            ("task_id", data.task_id),
            ("params", params),
            ("specimen_ids", specimen_ids),
            ("gene_id", data.gene_id),
            ("desc", "基因详情页-样本序列"),
            ("created_ts", datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        ]
        main_id = Dna(db).insert_main_table(collection="sg_gene_seq", data=mongo_data)
        Dna(db).update_db_record(collection="sg_gene_seq", query_dict={"_id": main_id},
                                 update_dict={"main_id": main_id})
        update_info = {str(main_id): "sg_gene_seq"}
        options = {
            "ref_fa": ref_fa,
            "ref_gff": ref_gff,
            "ssr_path": ssr_path,
            "specimen_ids": ";".join(specimen_ids),
            "pop_final_vcf": pop_final_vcf,
            "gene_id": data.gene_id,
            "update_info": json.dumps(update_info),
            "main_id": str(main_id),
            "is_wgs_result": is_wgs_result,
            "gene_name": data.gene_name
        }
        if hasattr(data, "project_type"):
            options["project_type"] = data.project_type
        main_table_name = "GeneDetail/SampleSeq_" + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3])
        self.set_sheet_data(name="dna_evolution.gene_sample_seq", member_id=member_id, project_sn=project_sn, task_id=data.task_id,
                            main_table_name=main_table_name, options=options, module_type="module", params=params, db_type=db)
        task_info = super(GeneSampleSeqAction, self).POST()
        task_info['id'] = str(main_id)
        return json.dumps(task_info)

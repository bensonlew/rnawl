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
        params = ["gene_id", "group_id", "group_dict", "task_id", "task_type", "submit_location"]
        for param in params:
            if not hasattr(data, param):
                var = []
                var.append(param)
                info = {"success": False, "info": "缺少%s参数!" % param, "code": "C3100601", "variables": var}
                return json.dumps(info)
        group_dict_ = json.loads(data.group_dict)
        specimen_ids = []
        for k in group_dict_.keys():
            for s in group_dict_[k]:
                specimen_ids.append(s)
        specimen_ids = list(set(specimen_ids))
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
                    "code": "C3100602", "variables": ""}
            return json.dumps(info)
        if pop_final_vcf.startswith("rerewrweset"):
            base_path = '/mnt/ilustre/data/' if str(data.client) == 'client01' else "/mnt/ilustre/tsanger-data/"
            pop_final_vcf_ = os.path.join(base_path, pop_final_vcf)
        elif pop_final_vcf.startswith("//"):
            region = result['region']  # 's3:'
            pop_final_vcf_ = ":".join([region.rstrip(':'), pop_final_vcf])
        elif pop_final_vcf.startswith("/mnt") or re.match(".*://.*", pop_final_vcf):
            pop_final_vcf_ = pop_final_vcf
        else:
            var = []
            var.append(pop_final_vcf)
            info = {"success": False, "info": "存入的文件%s路径格式不正确！" % pop_final_vcf,
                    "code": "C3100603", "variables": var}
            return json.dumps(info)
        ref_result = Dna("dna_wgs").find_one(collection="sg_species_version", query_dic={"_id": genome_version_id})
        # noinspection PyBroadException
        try:
            ref_fa = ref_result['ref']
            ref_gff = ref_result["gff"]
        except:
            info = {"success": False, "info": "sg_species_version表里没有ref、gff信息，请检查!",
                    "code": "C3100604", "variables": ""}
            return json.dumps(info)
        params_json = {
            "gene_id": data.gene_id,
            "group_id": data.group_id,
            "group_dict": group_dict_,
            "task_type": int(data.task_type),
            "submit_location": data.submit_location
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
            "specimen_ids": ";".join(specimen_ids),
            "pop_final_vcf": pop_final_vcf_,
            "gene_id": data.gene_id,
            "update_info": json.dumps(update_info),
            "main_id": str(main_id)
        }
        if hasattr(data, "project_type"):
            options["project_type"] = data.project_type
        main_table_name = "GeneDetail/SamplesSeq_" + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3])
        self.set_sheet_data(name="wgs.gene_samples_seq", member_id=member_id, project_sn=project_sn, task_id=data.task_id,
                            main_table_name=main_table_name, options=options, module_type="tool", params=params, db_type=db)
        task_info = super(GeneSampleSeqAction, self).POST()
        task_info['id'] = str(main_id)
        return json.dumps(task_info)

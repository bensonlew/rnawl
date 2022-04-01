# -*- coding: utf-8 -*-
# __author__ = 'zhaobinbin'
# modified 2018.01.11

import os
import re
import web
import json
import datetime
from mainapp.libs.signature import check_sig
from mainapp.models.mongo.dna import Dna
from mainapp.controllers.project.dna_controller import DnaController
from bson.objectid import ObjectId


class SnpCompareAction(DnaController):
    """
     变异位点比较分析接口
    """
    def __init__(self):
        super(SnpCompareAction, self).__init__(instant=False)

    @check_sig
    def  POST(self):
        data = web.input()
        print data
        params = ["task_id", "task_type", "analysis_model", "alle_number",
                  "submit_location", "chongmingming_result"]
        """
        analysis_model有两种情况，即single和multiple。
        
        """
        db = "dna_noref_wgs"
        result = Dna(db).find_one(collection="sg_task", query_dic={"task_id": data.task_id})
        try:
            vcf_path = result['pop_final_vcf']
            project_sn = result["project_sn"]
            member_id = result["member_id"]
            tag_file = result["populations_tag"]
        except:
            info = {"success": False, "info": "sg_task表里没有project_sn, member_id信息，请检查!", "code": "C3200205",
                    "variables": [""]}
            return json.dumps(info)
        vcf_path = Dna("dna_noref_wgs").set_file_path(data.task_id, vcf_path, data.client)
        tag_file = Dna("dna_noref_wgs").set_file_path(data.task_id, tag_file, data.client)
        if data.analysis_model == "single":  # 当为单条件的时候
            params.extend(["sample", "genotype", "marktype", "dep"])  # marktype即为diff或者是same
            for param in params:
                if not hasattr(data, param):
                    info = {"success": False, "info": "缺少%s参数!" , "variables": [param], "code" : "C3300211"}
                    return json.dumps(info)
            if data.marktype not in ["same", "diff"]:
                info = {"success": False, "info": "相同与不同的类型%s不合法!, 必须为same或者diff" , "variables":[data.marktype], "code" : "C3300212"}
                return json.dumps(info)
            if data.genotype not in ["homo", "hete", "all"]:
                info = {"success": False, "info": "基因型杂合纯合%s不合法!, 必须为homo或者heter或者all" , "variables": [data.genotype], "code" : "C3300213"}
                return json.dumps(info)
            params_json = {
                "sample": json.loads(data.sample),
                "marktype": data.marktype,
                "genotype": data.genotype,
                "dep": data.dep,
                "task_type": int(data.task_type),  # 用于表示是workflow，还是module还是tool之类的。
                "project_sn": project_sn,
                "submit_location":data.submit_location,
                "chongmingming_result": data.chongmingming_result,
                "analysis_model": data.analysis_model,
                "alle_number": data.alle_number
            }
            params = json.dumps(params_json, sort_keys=True, separators=(',', ':'))
            main_table_name = Dna("dna_noref_wgs").set_main_table_name("snp_compare", data.chongmingming_result)
            sub_name = []
            for sample in json.loads(data.sample):
                sample1 = sample.strip().split("|")[0]
                sample2 = sample.strip().split("|")[1]
                name = "{}_vs_{}_compare".format(sample1, sample2)
                sub_name.append(name)
            mongo_data = [
                ("subname", json.dumps(sub_name)),
                ("name", main_table_name),
                ("status", "start"),
                ("project_sn", project_sn),
                ("task_id", data.task_id),
                ("params", params),
                ("desc", "snp标记筛选！"),
                ("type", "single"),
                ("created_ts", datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
            ]
            collection = "sg_snp_compare"
            main_id = Dna(db).insert_main_table(collection=collection, data=mongo_data)
            Dna(db).update_db_record(collection=collection, query_dict={"_id": main_id},
                                     update_dict={"main_id": main_id})
            update_info = {str(main_id): collection}
            options = {
                "sample": data.sample,
                "genotype": data.genotype,
                "marktype": data.marktype,
                "vcf_file": vcf_path,
                "tag_file": tag_file,
                "dep": data.dep,
                "update_info": json.dumps(update_info),
                "main_id": str(main_id),
                "task_id": data.task_id,
                "analysis_model": data.analysis_model,
                "project_sn": project_sn,
                "alle_number": data.alle_number
            }
            self.set_sheet_data(name="noref_wgs.report.snp_compare", member_id=member_id, project_sn=project_sn,
                                task_id=data.task_id, main_table_name="snp_compare/" + main_table_name,
                                options=options, module_type="workflow",   params=params, db_type=db)
            task_info = super(SnpCompareAction, self).POST()
            task_info['id'] = str(main_id)
            return json.dumps(task_info)

        elif data.analysis_model == "multiple":  # 当为多条件的时候
            sample = json.loads(data.sample)
            if sample == [""]:
                pass
            else:
                for i in json.loads(data.marktype):
                    if i not in ["diff", "same"]:
                        info = {"success": False, "info": "相同与不同的类型%s不合法!, 必须为diff或者same" , "variables": [data.marktype], "code" : "C3300214"}
                        return json.dumps(info)
            group_id = json.loads(data.group_id)
            group_list = []
            if group_id == [""]:
                pass
            else:
                for id in group_id:
                    if id == "":
                        info = {"success": False, "info": "%s存在group_id为空的情况，请核查!" , "variables": [data.group_id], "code" : "C3300215"}
                        return json.dumps(info)
                group_dict_new = json.loads(data.group_dict)
                for dict1 in group_dict_new:
                    for key in dict1.keys():
                        group_list.append(key + ":" + ",".join(dict1[key]))
            params_json = {
                "group_ids": data.group_ids,
                "cate_groups": data.cate_groups,
                "spe_groups": data.spe_groups,
                "group_name": data.group_name,
                "group_id": json.loads(data.group_id),
                "group_dict": json.loads(data.group_dict),
                "maf": json.loads(data.maf),
                "ad": json.loads(data.ad),
                "max_miss": json.loads(data.max_miss),
                "sample": json.loads(data.sample),
                "marktype": json.loads(data.marktype),
                "genotype": json.loads(data.genotype),
                "dep": json.loads(data.dep),
                "task_type": int(data.task_type),  # 用于表示是workflow，还是module还是tool之类的。
                "analysis_model": data.analysis_model,
                "project_sn": project_sn,
                "submit_location": data.submit_location,
                "chongmingming_result": data.chongmingming_result,
                "alle_number": data.alle_number
            }
            params = json.dumps(params_json, sort_keys=True, separators=(',', ':'))
            main_table_name = Dna("sg_snp_compare").set_main_table_name("snp_compare", data.chongmingming_result)
            mongo_data = [
                ("name", main_table_name),
                ("status", "start"),
                ("project_sn", project_sn),
                ("task_id", data.task_id),
                ("params", params),
                ("type", "multiple"),
                ("desc", "multiple"),
                ("created_ts", datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
            ]
            main_id = Dna(db).insert_main_table(collection="sg_snp_compare", data=mongo_data)
            Dna(db).update_db_record(collection="sg_snp_compare", query_dict={"_id": main_id},
                                     update_dict={"main_id": main_id})
            update_info = {str(main_id): "sg_snp_compare"}
            options = {
                "group": json.dumps(group_list),
                "analysis_model": data.analysis_model,
                "maf": data.maf,
                "ad": data.ad,
                "max_miss": data.max_miss,
                "update_info": json.dumps(update_info),
                "main_id": str(main_id),
                "task_id": data.task_id,
                "vcf_file": vcf_path,
                "tag_file": tag_file,
                "sample": data.sample,
                "marktype": data.marktype,
                "dep": data.dep,
                "genotype": data.genotype,
                "project_sn": project_sn,
                "alle_number": data.alle_number
            }
            self.set_sheet_data(name="noref_wgs.report.snp_compare", member_id=member_id, project_sn=project_sn,
                                task_id=data.task_id,
                                main_table_name=main_table_name, options=options, module_type="workflow", params=params,
                                db_type=db)
            task_info = super(SnpCompareAction, self).POST()
            task_info['id'] = str(main_id)
            return json.dumps(task_info)

# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# modified 2018.04.14

import os
import re
import web
import json
import datetime
from mainapp.libs.signature import check_sig
from mainapp.models.mongo.dna import Dna
from mainapp.controllers.project.dna_controller import DnaController


class DoubleGroupCompareAction(DnaController):
    """
    snp/indel样本组与样本组的比较分析接口
    """
    def __init__(self):
        super(DoubleGroupCompareAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        print data
        required_params = ["group_id", "group_dict", "compare_group", "location", "funtype", "efftype", "maf1", "maf2", "dep1",
                           "dep2", "miss1", "miss2", "len1", "len2", "task_id", "task_type", "submit_location"]
        for param in required_params:
            if not hasattr(data, param):
                var = []
                var.append(param)
                info = {"success": False, "info": "缺少%s参数!" % param, "code": "C3100401", "variables": var}
                return json.dumps(info)
        group_dict_ = {}
        try:
            group1_name = data.compare_group.split("|")[0]
            group2_name = data.compare_group.split("|")[1]
        except:
            var = []
            var.append(data.compare_group)
            info = {"success": False, "info": "%s差异分组传入有误，请检查" % data.compare_group,
                    "code": "C3100402", "variables": var}
            return json.dumps(info)
        try:
            group_dict_ = json.loads(data.group_dict)
            group1_samples = group_dict_[group1_name]
            group2_samples = group_dict_[group2_name]
            group1 = group1_name + ":" + ','.join(group1_samples)
            group2 = group2_name + ":" + ','.join(group2_samples)
        except:
            var = []
            var.append(data.group_dict)
            info = {"success": False, "info": "%s分组传入有误，请检查" % data.group_dict,
                    "code": "C3100403", "variables": var
                    }
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
        except:
            info = {"success": False, "info": "sg_task表里没有project_sn, member_id,pop_final_vcf信息，请检查!",
                    "code": "C3100404", "variables": ""}
            return json.dumps(info)
        if pop_final_vcf.startswith("rerewrweset"):
            base_path = '/mnt/ilustre/data/' if str(data.client) == 'client01' else "/mnt/ilustre/tsanger-data/"
            pop_final_vcf_ = os.path.join(base_path, pop_final_vcf)
        elif pop_final_vcf.startswith("//"):
            region = result['region']   # 's3:'
            pop_final_vcf_ = ":".join([region.rstrip(':'), pop_final_vcf])
        elif pop_final_vcf.startswith("/mnt") or re.match(".*://.*", pop_final_vcf):
            pop_final_vcf_ = pop_final_vcf
        else:
            var = []
            var.append(pop_final_vcf)
            info = {"success": False, "info": "存入的文件%s路径格式不正确！" % pop_final_vcf,
                    "code": "C3100405", "variables": ""}
            return json.dumps(info)
        # if not os.path.exists(pop_final_vcf):
        #     info = {"success": False, "info": "sg_task表里的pop_final_vcf文件不存在，请检查!"}
        #     return json.dumps(info)
        if data.location == ",," or data.location == ",":
            location = None
        else:
            location = data.location
        if re.match(r".+,,", data.location):
            info = {"success": False, "info": "基因组区域不能只填染色体编号!", "code": "C3100406", "variables": ""}
            return json.dumps(info)
        funtype = data.funtype if data.funtype != "," else None
        efftype = data.efftype if data.efftype != "," else None
        maf1 = data.maf1 if data.maf1 != "," else None
        maf2 = data.maf2 if data.maf2 != "," else None
        dep1 = data.dep1 if data.dep1 != "," else None
        dep2 = data.dep2 if data.dep2 != "," else None
        miss1 = data.miss1 if data.miss1 != "," else None
        miss2 = data.miss2 if data.miss2 != "," else None
        len1 = data.len1 if data.len1 else "-1"
        len2 = data.len2 if data.len2 else "100000000"
        len = len1 + "," + len2
        print len
        if data.submit_location == "snpcompare_sgroup":
            analysis_type = "snp"
            if data.len1 != "1" or data.len2 != "1":
                info = {"success": False, "info": "为snp分析时len1和len2需为1!", "code": "C3100407", "variables": ""}
                return json.dumps(info)
        else:
            analysis_type = "indel"
        params_json = {
            "group_id": data.group_id,
            "group_dict": group_dict_,
            "compare_group": data.compare_group,
            "location": data.location,
            "funtype": data.funtype,
            "efftype": data.efftype,
            "maf1": data.maf1,
            "maf2": data.maf2,
            "dep1": data.dep1,
            "dep2": data.dep2,
            "miss1": data.miss1,
            "miss2": data.miss2,
            "len1": data.len1,
            "len2": data.len2,
            "task_type": int(data.task_type),
            "submit_location": data.submit_location
        }
        params = json.dumps(params_json, sort_keys=True, separators=(',', ':'))
        mongo_data = [
            ("name", "{}_vs_{}_{}_".format(group1_name, group2_name, analysis_type) +
             str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S"))),
            ("status", "start"),
            ("project_sn", project_sn),
            ("task_id", data.task_id),
            ("params", params),
            ("type", "two_group"),
            ("desc", "样本组与样本组差异比较分析"),
            ("created_ts", datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        ]
        if analysis_type == "snp":
            main_table_name = "sg_snp_compare"
            main_id = Dna(db).insert_main_table(collection="sg_snp_compare", data=mongo_data)
            Dna(db).update_db_record(collection="sg_snp_compare", query_dict={"_id": main_id},
                                            update_dict={"main_id": main_id})
            update_info = {str(main_id): "sg_snp_compare"}
        else:
            main_table_name = "sg_indel_compare"
            main_id = Dna(db).insert_main_table(collection="sg_indel_compare", data=mongo_data)
            Dna(db).update_db_record(collection="sg_indel_compare", query_dict={"_id": main_id},
                                            update_dict={"main_id": main_id})
            update_info = {str(main_id): "sg_indel_compare"}
        options = {
            "pop_final_vcf": pop_final_vcf_,
            "group1": group1,
            "group2": group2,
            "analysis_type": analysis_type,
            "location": location,
            "funtype": funtype,
            "efftype": efftype,
            "maf1": maf1,
            "maf2": maf2,
            "dep1": dep1,
            "dep2": dep2,
            "miss1": miss1,
            "miss2": miss2,
            "len1": len,
            "len2": len,
            "update_info": json.dumps(update_info),
            "main_id": str(main_id),
            "task_id": data.task_id
        }
        if hasattr(data, "project_type"):
            options["project_type"] = data.project_type
        if analysis_type == "snp":
            main_table_name = "SNPCompare/{}_vs_{}_{}_".format(group1_name, group2_name, analysis_type) + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S"))
        else:
            main_table_name = "IndelCompare/{}_vs_{}_{}_".format(group1_name, group2_name, analysis_type) + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S"))
        self.set_sheet_data(name="wgs.report.double_group_compare", member_id=member_id, project_sn=project_sn, task_id=data.task_id,
                            main_table_name=main_table_name, options=options, module_type="workflow", params=params, db_type=db)
        task_info = super(DoubleGroupCompareAction, self).POST()
        task_info['id'] = str(main_id)
        return json.dumps(task_info)

# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# modified 2018.04.16

import os
import re
import web
import json
import datetime
from mainapp.libs.signature import check_sig
from mainapp.models.mongo.dna import Dna
from mainapp.controllers.project.dna_controller import DnaController


class SingleGroupCompareAction(DnaController):
    """
    snp/indel组内样本间的比较分析接口
    """
    def __init__(self):
        super(SingleGroupCompareAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        print data
        required_params = ["group_id", "group_dict", "location", "funtype", "efftype", "maf", "dep", "miss",
                           "len1", "len2", "task_id", "task_type", "submit_location"]
        for param in required_params:
            if not hasattr(data, param):
                var = []
                var.append(param)
                info = {"success": False, "info": "缺少%s参数!" % param, "code": "C3101101", "variables": var}
                return json.dumps(info)
        # noinspection PyBroadException
        try:
            group_dict_ = json.loads(data.group_dict)
            group_name = group_dict_.keys()[0]
            group_samples = group_dict_[group_name]
            group = group_name + ":" + ','.join(group_samples)
        except:
            var = []
            var.append(data.group_dict)
            info = {"success": False, "info": "{}分组传入有误，请检查" % data.group_dict, "code": "C3101102", "variables": var}
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
                    "code": "C3101103", "variables": ""}
            return json.dumps(info)
        pop_final_vcf = Dna(db).set_file_path(data.task_id, pop_final_vcf, data.client)
        # if not os.path.exists(pop_final_vcf):
        #     info = {"success": False, "info": "sg_task表里的pop_final_vcf文件不存在，请检查!"}
        #     return json.dumps(info)
        if re.match(r".+,,", data.location):
            info = {"success": False, "info": "基因组区域不能只填染色体编号!",
                    "code": "C3101105", "variables": ""}
            return json.dumps(info)
        if data.location == ",," or data.location == ",":
            location = None
        else:
            location = data.location
        s_funtype = []
        for s in data.funtype.split(","):
            if s:
                s_funtype.append(s)
        funtype = ','.join(s_funtype) if s_funtype else None
        s_efftype = []
        for s in data.efftype.split(","):
            if s:
                s_efftype.append(s)
        efftype = ','.join(s_efftype) if s_efftype else None
        maf = data.maf if data.maf != "," else None
        dep = data.dep if data.dep != "," else None
        miss = data.miss if data.miss != "," else None
        len1 = data.len1 if data.len1 else "-1"
        len2 = data.len2 if data.len2 else "100000000"
        len = len1 + "," + len2
        print len
        if data.submit_location == "snpcompare_group":
            analysis_type = "snp"
            if data.len1 != "1" or data.len2 != "1":
                info = {"success": False, "info": "为snp分析时len1和len2需为1!",
                        "code": "C3101106", "variables": ""}
                return json.dumps(info)
        else:
            analysis_type = "indel"
        params_json = {
            "group_id": data.group_id,
            "group_dict": group_dict_,
            "location": data.location,
            "funtype": data.funtype,
            "efftype": data.efftype,
            "maf": data.maf,
            "dep": data.dep,
            "miss": data.miss,
            "len1": data.len1,
            "len2": data.len2,
            "task_type": int(data.task_type),
            "submit_location": data.submit_location
        }
        params = json.dumps(params_json, sort_keys=True, separators=(',', ':'))
        mongo_data = [
            ("name", "{}_{}_".format(group_name, analysis_type) +
             str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f"))),
            ("status", "start"),
            ("project_sn", project_sn),
            ("task_id", data.task_id),
            ("params", params),
            ("type", "one_group"),
            ("desc", "组内样本间差异比较分析"),
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
            "pop_final_vcf": pop_final_vcf,
            "group": group,
            "analysis_type": analysis_type,
            "location": location,
            "funtype": funtype,
            "efftype": efftype,
            "maf": maf,
            "dep": dep,
            "miss": miss,
            "len": len,
            "update_info": json.dumps(update_info),
            "main_id": str(main_id),
            "task_id": data.task_id
        }
        if hasattr(data, "project_type"):
            options["project_type"] = data.project_type
        if analysis_type == "snp":
            main_table_name = "SNPCompare/{}_{}_".format(group_name, analysis_type) + \
                              str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f"))
        else:
            main_table_name = "IndelCompare/{}_{}_".format(group_name, analysis_type) + \
                              str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f"))
        self.set_sheet_data(name="wgs.report.single_group_compare", member_id=member_id, project_sn=project_sn,
                            task_id=data.task_id, main_table_name=main_table_name, options=options,
                            module_type="workflow", params=params, db_type=db)
        task_info = super(SingleGroupCompareAction, self).POST()
        task_info['id'] = str(main_id)
        return json.dumps(task_info)

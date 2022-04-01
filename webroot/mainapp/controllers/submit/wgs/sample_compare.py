# -*- coding: utf-8 -*-
# __author__ = 'HONGDONG'
# modified 20180414

import os
import re
import web
import json
import datetime
from mainapp.libs.signature import check_sig
from mainapp.models.mongo.dna import Dna
from mainapp.controllers.project.dna_controller import DnaController


class SampleCompareAction(DnaController):
    """
     样本比较分析接口， 当该接口用于snp的时候len1与len2默认都为1
    """
    def __init__(self):
        super(SampleCompareAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        print data
        params = ["sample1", "sample2", "is_same", "funtype", "efftype", "len1", "len2", "location",
                  "task_id", "task_type", "submit_location"]
        for param in params:
            if not hasattr(data, param):
                var = []
                var.append(param)
                info = {"success": False, "info": "缺少%s参数!" % param, "code": "C3100901", "variables": var}
                return json.dumps(info)
        if data.is_same not in ["true", "false"]:
            var = []
            var.append(data.is_same)
            info = {"success": False, "info": "相同与不同的类型%s不合法!, 必须为true或者false" % data.is_same,
                    "code": "C3100902", "variables": var}
            return json.dumps(info)
        if data.submit_location not in ["snpcompare_specimen", "indelcompare_specimen"]:
            var = []
            var.append(data.submit_location)
            info = {"success": False, "info": "submit_location的类型%s不合法!, 因为我们根据该字段进行区分snp与indel分析，"
                                              "必须为snpcompare_specimen或者indelcompare_specimen" % data.submit_location,
                    "code": "C3100903", "variables": var}
            return json.dumps(info)
        if "dep1" in data.keys():
            dep1 = data.dep1
        else:
            dep1 = ","
        if "dep2" in data.keys():
            dep2 = data.dep2
        else:
            dep2 = ","
        if re.match(r".+,,", data.location):  # 增加location的限制，modified by zengjing 20180612
            info = {"success": False, "info": "基因组区域不能只填染色体编号!", "code": "C3100904", "variables": ""}
            return json.dumps(info)
        params_json = {
            "sample1": data.sample1,
            "sample2": data.sample2,
            "is_same": data.is_same,
            "funtype": data.funtype,
            "efftype": data.efftype,
            "len1": data.len1,
            "len2": data.len2,
            "dep1": dep1,
            "dep2": dep2,
            "location": data.location,
            "task_type": int(data.task_type),
            "submit_location": data.submit_location
        }
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
            info = {"success": False, "info": "sg_task表里没有project_sn, member_id信息，请检查!",
                    "code": "C3100905", "variables": ""}
            return json.dumps(info)
        pop_final_vcf = Dna(db).set_file_path(data.task_id, pop_final_vcf, data.client)
        # if not os.path.exists(pop_final_vcf):
        #     info = {"success": False, "info": "{}文件不存在！".format(pop_final_vcf)}
        #     return json.dumps(info)
        params = json.dumps(params_json, sort_keys=True, separators=(',', ':'))
        if data.submit_location == "indelcompare_specimen":  # 完善indel与snp的区分
            analysis_type = 'indel'
            if data.len1 and int(data.len1) > 1:
                len1 = data.len1
            else:
                len1 = '2'
        else:
            analysis_type = 'snp'
            len1 = '1'
        main_table_name = '{}_{}_{}_compare_'.format(data.sample1, data.sample2, analysis_type) + \
                          datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        mongo_data = [
            ("name", main_table_name),
            ("status", "start"),
            ("project_sn", project_sn),
            ("task_id", data.task_id),
            ("params", params),
            ("desc", "样本比较分析！"),
            ("type", "sample"),
            ("created_ts", datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        ]
        collection = "sg_indel_compare" if data.submit_location == "indelcompare_specimen" else "sg_snp_compare"
        main_id = Dna(db).insert_main_table(collection=collection, data=mongo_data)
        Dna(db).update_db_record(collection=collection, query_dict={"_id": main_id},
                                 update_dict={"main_id": main_id})
        update_info = {str(main_id): collection}
        options = {
            "sample1": data.sample1,
            "sample2": data.sample2,
            "is_same": data.is_same,
            "vcf_file": pop_final_vcf,
            "funtype": data.funtype,
            "efftype": data.efftype,
            "len1": len1,
            "len2": data.len2,
            "dep1": dep1,
            "dep2": dep2,
            "location": data.location,
            "update_info": json.dumps(update_info),
            "main_id": str(main_id),
            "submit_location": data.submit_location,
            "task_id": data.task_id
        }
        if hasattr(data, "project_type"):
            options["project_type"] = data.project_type
        self.set_sheet_data(name="wgs.report.sample_compare", member_id=member_id, project_sn=project_sn,
                            task_id=data.task_id, main_table_name="sample_compare/" + main_table_name,
                            options=options, params=params, db_type=db)
        task_info = super(SampleCompareAction, self).POST()
        task_info['id'] = str(main_id)
        return json.dumps(task_info)

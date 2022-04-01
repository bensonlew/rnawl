# -*- coding: utf-8 -*-
# __author__ = 'Binbin Zhao
# modified 20180904


import web
import json
import datetime
import os
from bson.objectid import ObjectId
from types import StringTypes
from mainapp.libs.signature import check_sig
from mainapp.models.mongo.dna import Dna
from mainapp.controllers.project.dna_controller import DnaController


class LdAnalysisAction(DnaController):
    """
    群体进化 基因注释接口
    """

    def __init__(self):
        super(LdAnalysisAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        print data
        params = [
            "task_id",
            "task_type",
            "group_id",
            "group_dict",
            "max_missing",
            "min_maf",
            "max_maf",
           "vcf_id", "submit_location", "chongmingming_result", "maxdp", "mindp"]
        for param in params:
            if not hasattr(data, param):
                var = []
                var.append(param)
                info = {"success": False, "info": "缺少%s参数!" % param, "code": "C3200701", "variables": var}
                return json.dumps(info)
        dict1 = json.loads(data.group_dict)
        new_dict = {}
        for value in dict1.keys():
            new_dict[value] = ",".join(dict1[value])
        if data.max_missing == "":
            data.max_missing = None
        if data.min_maf == "":
            data.min_maf = None
        # try:
        #     min_maf = data.maf.split(",")[0]
        #     max_maf = data.maf.split(",")[1]
        #     if min_maf == "":
        #         min_maf = None
        #     if max_maf == "":
        #         max_maf = None
        # except:
        #     info = {"success": False, "info": "maf格式不正确", "code": "C3200702", "variables": ""}
        #     return json.dumps(info)
        # try:
        #     min_dp = data.average_dp.split(",")[0]
        #     max_dp = data.average_dp.split(",")[1]
        #     if min_dp == "":
        #         min_dp = None
        #     if max_dp == "":
        #         max_dp = None
        # except:
        #     info = {"success": False, "info": "average_dp格式不正确", "code": "C3200703", "variables": ""}
        #     return json.dumps(info)
        if data.maxdp == "":
            data.maxdp = None
        if data.mindp == "":
            data.mindp = None
        task_result = Dna("dna_evolution").find_one(collection="sg_task", query_dic={"task_id": data.task_id})
        if not task_result:
            var = []
            var.append(data.task_id)
            info = {"success": False,
                    "info": "sg_task表里没有task_id: %s对应的信息，请检查!" % data.task_id, "code": "C3200704",
                    "variables": var}
            return json.dumps(info)
        project_sn = task_result["project_sn"]
        member_id = task_result["member_id"]
        vcf_file = Dna("dna_evolution").find_one(collection="sg_variant_compare_filter",
                                                 query_dic={"_id": ObjectId(data.vcf_id)})
        if not vcf_file:
            var = []
            var.append(data.vcf_id)
            info = {"success": False,
                    "info": "sg_variant_compare_filter表里没有_id: %s对应的信息，请检查!" % data.vcf_id,
                    "code": "C3200705", "variables": var}
            return json.dumps(info)
        else:
            # noinspection PyBroadException
            try:
                vcf_file_path = vcf_file['vcf_path']
            except:
                info = {"success": False,
                        "info": "sg_variant_compare_filter中没有vcf_path字段，请核查！", "code": "C3200706",
                        "variables": ""}
                return json.dumps(info)
        vcf_file_path = Dna("dna_evolution").set_file_path(data.task_id, vcf_file_path, data.client)
        # project_sn = "test_ld_decay"  # 测试的时候使用
        # member_id = "12dfsfdsfdsfksfls"  # 测试的时候使用
        # file_path = "/mnt/ilustre/users/sanger-dev/sg-users/zhaobinbin/GeneticEvolution/ld_decay" \
        #             "/pop1.sort.final.vcf"  # 测试专用路径
        # file_path = "/mnt/ilustre/users/sanger-dev/sg-users/xuanhongdong/evolution/test_data/vcf_data_201809214/" \
        #             "1.vcf"
        # main_id = "5a8e7a56e95d0d5f827bd7fe"
        params_json = {
            "task_type": int(data.task_type),
            "group_id": data.group_id,
            "group_dict": json.loads(data.group_dict),
            "max_missing": data.max_missing,
            "min_maf": data.min_maf,
            "max_maf": data.max_maf,
            "min_dp": data.mindp,
            "max_dp": data.maxdp,
            "submit_location": data.submit_location,
            "chongmingming_result": data.chongmingming_result
        }
        params = json.dumps(
            params_json,
            sort_keys=True,
            separators=(
                ',',
                ':'))     # 返回json
        # main_table_name = 'Ld_Analysis_' + \
        #     datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")
        main_table_name = Dna("dna_evolution").set_main_table_name("LdAnalysis", data.chongmingming_result)
        mongo_data = [
            ("name", main_table_name),
            ("project_sn", project_sn),
            ("task_id", data.task_id),
            ("created_ts", datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ("status", "start"),
            ("params", params),
            ("desc", "连锁不平衡"),
        ]
        main_id = Dna("dna_evolution").insert_main_table(
            collection="sg_ld", data=mongo_data)
        Dna("dna_evolution").update_db_record(
            collection="sg_ld", query_dict={
                "_id": main_id}, update_dict={
                "main_id": main_id})
        update_info = {str(main_id): "sg_ld"}
        options = {
            "task_id": data.task_id,  # 页面另外两个信息前端可自查表，不用存params?!
            "group_dict": json.dumps(new_dict),
            "max_missing": data.max_missing,
            "min_maf": data.min_maf,
            "max_maf": data.max_maf,
            "min_dp": data.mindp,
            "max_dp": data.maxdp,
            "project_sn": project_sn,
            "main_id": str(main_id),
            "update_info": json.dumps(update_info),     # 给work_flow用
            "vcf_file": vcf_file_path
        }
        self.set_sheet_data(name="dna_evolution.ld_analysis", member_id=member_id, project_sn=project_sn,
                            task_id=data.task_id, main_table_name="ld/" + main_table_name,
                            options=options, module_type="module", params=params,
                            db_type="dna_evolution", analysis_name="LdAnalysis")
        task_info = super(LdAnalysisAction, self).POST()
        task_info['id'] = str(main_id)
        return json.dumps(task_info)

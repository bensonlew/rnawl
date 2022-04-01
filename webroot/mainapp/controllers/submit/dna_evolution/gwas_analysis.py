# -*- coding: utf-8 -*-
# __author__ = 'qing_mei
# modified 20180825
# controller.submit

import web
import json
import datetime
import os
from bson.objectid import ObjectId
from types import StringTypes
from mainapp.libs.signature import check_sig
from mainapp.models.mongo.dna import Dna
from mainapp.controllers.project.dna_controller import DnaController


class GwasAnalysisAction(DnaController):
    """
    群体进化：GWAS关联分子接口
    """
    def __init__(self):
        super(GwasAnalysisAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        print('开始GWAS关联分析')
        print data
        params = ["vcf_id", "upload_trit_id", "max_missing", "min_maf", "max_maf", "min_dp", "max_dp",
                  "chongmingming_result", "submit_location", "task_id"]
        for param in params:
            if not hasattr(data, param):
                info = {"success": False, "info": "缺少%s参数!" % param}
                return json.dumps(info)
        if data.upload_trit_id == "":
            info = {"success": False, "info": "%s参数为空，请核查!" % data.upload_trit_id}
            return json.dumps(info)
        if data.chrs_list == "":
            info = {"success": False, "info": "%s参数为空，请核查!" % data.chrs_list}
            return json.dumps(info)
        upload_trait_result = Dna("dna_evolution").find_one(collection="sg_upload_trait",
                                                            query_dic={"_id": ObjectId(data.upload_trit_id)})
        upload_trait_path = ""
        if "file_url" in upload_trait_result.keys():
            upload_trait_path = upload_trait_result["file_url"]
        else:
            info = {"success": False, "info": "sg_upload_trait表缺少trait_path：%s参数!" % upload_trait_path}
            return json.dumps(info)
        # upload_trait_path = self.get_target_path(data.client, upload_trait_path)
        task_result = Dna("dna_evolution").find_one(collection="sg_task", query_dic={"task_id": data.task_id})
        if not task_result:
            info = {"success": False,
                    "info": "sg_task表里没有task_id: %s对应的信息，请检查!" % data.task_id}
            return json.dumps(info)
        project_sn = task_result["project_sn"]
        member_id = task_result["member_id"]
        vcf_file = Dna("dna_evolution").find_one(collection="sg_variant_compare_filter",
                                                 query_dic={"_id": ObjectId(data.vcf_id)})
        if not vcf_file:
            info = {"success": False,
                    "info": "sg_variant_compare_filter表里没有_id: %s对应的信息，请检查!" % data.vcf_id}
            return json.dumps(info)
        else:
            # noinspection PyBroadException
            try:
                vcf_file_path = vcf_file['vcf_path']
            except:
                info = {"success": False,
                        "info": "sg_variant_compare_filter中没有vcf_path字段，请核查！"}
                return json.dumps(info)
        vcf_file_path = Dna("dna_evolution").set_file_path(data.task_id, vcf_file_path, data.client)
        params_json = {
            "vcf_id": data.vcf_id,
            "upload_trit_id": data.upload_trit_id,
            "max_missing": data.max_missing,
            "min_maf": data.min_maf,
            "max_maf": data.max_maf,
            "min_dp": data.min_dp,
            "max_dp": data.max_dp,
            "task_id": data.task_id,
            "submit_location": data.submit_location,
            "chrs_list": data.chrs_list
        }
        options = {
            "vcf_path": vcf_file_path,
            "max_missing": data.max_missing,
            "min_maf": data.min_maf,
            "max_maf": data.max_maf,
        }
        if data.min_dp:
            options.update({"min_dp": data.min_dp})
        if data.max_dp:
            options.update({"max_dp": data.max_dp})
        # detail_info_path = self.get_target_path(data.client, detail_info_path)
        # main_table_name = data.name if data.name != "" else 'GwasAnalysis_' + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")
        main_table_name = Dna("dna_evolution").set_main_table_name("GwasAnalysis", data.chongmingming_result)
        params = json.dumps(params_json, sort_keys=True, separators=(',', ':'))  # 返回json
        mongo_data = [
            ("name", main_table_name),
            ("project_sn", project_sn),
            ("task_id", data.task_id),
            ("created_ts", datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ("status", "start"),
            ("params", params),
            ("desc", "Gwas接口分析")
        ]
        main_id = Dna("dna_evolution").insert_main_table(collection="sg_gwas_analysis", data=mongo_data)
        Dna("dna_evolution").update_db_record(collection="sg_gwas_analysis", query_dict={"_id": main_id},
                                              update_dict={"main_id": main_id})
        update_info = {str(main_id): "sg_gwas_analysis"}
        result = Dna("dna_evolution").find_one(collection="sg_gwas_analysis", query_dic={"_id": main_id})
        if not result:
            info = {"success": False,
                    "info": "sg_gwas_analysis表里没有_id: %s对应的信息，请检查!" % main_id}
            return json.dumps(info)
        options.update({
            "main_id": str(main_id),
            "task_id": data.task_id,
            "vcf_path": vcf_file_path,    # ????
            "upload_trait_path": upload_trait_path,
            "update_info": json.dumps(update_info),
            "chrs_list": data.chrs_list
        })
        self.set_sheet_data(name="dna_evolution.gwas_analysis", member_id=member_id, project_sn=project_sn,
                            task_id=data.task_id, main_table_name="Gwas_analysis/" + main_table_name, options=options,
                            module_type="module", params=params, db_type="dna_evolution", analysis_name="Gwas")
        task_info = super(GwasAnalysisAction, self).POST()
        task_info['id'] = str(main_id)
        return json.dumps(task_info)

    def get_target_path(self, client, path):
        """
        获取远程磁盘的路径
        :return:
        """
        if client not in ['client01', 'client03']:
            raise Exception("client{}类型不正确！".format(client))
        if client == 'client01':
            target_path = os.path.join("/mnt/ilustre/data", path)
        else:
            target_path = os.path.join("/mnt/ilustre/tsanger-data", path)
        return target_path

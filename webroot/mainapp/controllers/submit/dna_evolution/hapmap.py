# -*- coding: utf-8 -*-
# __author__ = 'qing_mei
# modified 20180830
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


class HapmapAction(DnaController):
    """
    群体进化 单倍体图谱
    """

    def __init__(self):
        super(HapmapAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        print data
        params = ["sg_gwas_analysis_id", "chongmingming_result", "submit_location", "task_id"]
        for param in params:
            if not hasattr(data, param):
                info = {"success": False, "info": "缺少%s参数,请检查" % param}
                return json.dumps(info)
        if not hasattr(data, "p_value") and not hasattr(data, "q_value"):
            info = {"success": False, "info": "参数%s和%s至少要有一个，请检查" % ("p_value", "q_value")}
            return json.dumps(info)
        if not hasattr(data, "distance") and not hasattr(data, "ld_r2"):
            info = {"success": False, "info": "参数%s和%s至少要有一个，请检查" % ("distance", "ld_r2")}
            return json.dumps(info)
        if hasattr(data, "p_value"):
            if float(data.p_value) <= 0 or float(data.p_value) >= 1:
                info = {"success": False, "info": "P值:%s必需为(0,1),请检查" % data.p_value}
                return json.dumps(info)
        if hasattr(data, "ld_r2"):
            if float(data.ld_r2) <= 0 or float(data.ld_r2) >= 1:
                info = {"success": False, "info": "R2值:%s必需为(0,1),请检查" % data.ld_r2}
                return json.dumps(info)
        task_result = Dna("dna_evolution").find_one(collection="sg_task", query_dic={"task_id": data.task_id})
        if not task_result:
            info = {"success": False,
                    "info": "sg_task表里没有task_id: %s对应的信息，请检查!" % data.task_id}
            return json.dumps(info)
        project_sn = task_result["project_sn"]
        member_id = task_result["member_id"]
        # anno_summary_path = "/mnt/ilustre/users/sanger-dev/sg-users/liuwentian/evoluation/test_data_1023/anno.summary"
        genome_version_id = task_result["genome_version_id"]
        is_wgs_result = task_result["is_wgs_result"]
        mongo_dna = "dna_evolution"
        if is_wgs_result == "yes":
            mongo_dna = "dna_wgs"
        anno_summary = Dna(mongo_dna).find_one(collection="sg_species_version",
                                               query_dic={"_id": ObjectId(genome_version_id)})
        try:
            anno_summary_path = anno_summary['anno']
        except:
            info = {"success": False,
                    "info": "sg_species_version中没有anno字段，请核查！"}
            return json.dumps(info)
        # anno_summary_path = Dna("dna_evolution").set_file_path(data.task_id, anno_summary_path, data.client)
        # 此处需要根据task_id来查基因组数据库的表，然后得到anno.summary的文件rere路径，然后再用下面代码换成绝对路径，传给接口的程序
        # anno_summary_path =  self.get_target_path(data.client, anno_summary_path)
        gwas_result = Dna("dna_evolution").find_one(collection="sg_gwas_analysis",
                                                    query_dic={"_id": ObjectId(data.sg_gwas_analysis_id)})
        if not gwas_result:
            info = {"success": False,
                    "info": "sg_gwas_analysis表里没有_id: %s对应的信息，请检查!" % data.sg_gwas_analysis_id}
            return json.dumps(info)
        for i in ['csv_dir', 'recode_vcf_path', 'trait_list']:
            if i not in gwas_result.keys():
                info = {"success": False, "info": "sg_gwas_analysis表里没有{}，请检查!".format(i)}
                return json.dumps(info)
        csv_dir = gwas_result["csv_dir"]
        recode_vcf_path = gwas_result["recode_vcf_path"]
        trait_list = gwas_result["trait_list"]
        # csv_dir = self.get_target_path(data.client, csv_dir)
        # recode_vcf_path = self.get_target_path(data.client, recode_vcf_path)
        options = {}
        params_json = {
            "sg_gwas_analysis_id": data.sg_gwas_analysis_id,
            "submit_location": data.submit_location
        }
        if hasattr(data, "distance"):
            params_json['distance'] = int(data.distance) * 1000
            options['distance'] = int(data.distance) * 1000
        elif hasattr(data, "ld_r2"):
            params_json['ld_r2'] = float(data.ld_r2)
            options['ld_r2'] = float(data.ld_r2)
        if hasattr(data, "p_value"):
            params_json['p_value'] = float(data.p_value)
            options['p_value'] = float(data.p_value)
        elif hasattr(data, "q_value"):
            params_json['q_value'] = float(data.q_value)
            options['q_value'] = float(data.q_value)
        params = json.dumps(params_json, sort_keys=True, separators=(',', ':'))     # 返回json
        # main_table_name = 'Hapmap_' + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")
        main_table_name = Dna("dna_evolution").set_main_table_name("Hapmap", data.chongmingming_result)
        mongo_data = [
            ("name", main_table_name),
            ("project_sn", project_sn),
            ("task_id", data.task_id),
            ("created_ts", datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ("status", "start"),
            ("params", params),
            ("desc", "单倍体图谱接口"),
            ("trait_list", trait_list),
        ]
        main_id = Dna("dna_evolution").insert_main_table(collection="sg_hapmap", data=mongo_data)
        Dna("dna_evolution").update_db_record(collection="sg_hapmap", query_dict={"_id": main_id},
                                              update_dict={"main_id": main_id})
        update_info = {str(main_id): "sg_hapmap"}
        result = Dna("dna_evolution").find_one(collection="sg_hapmap", query_dic={"_id": main_id})
        if not result:
            info = {"success": False,
                    "info": "sg_hapmap表里没有main_id: %s对应的信息，请检查!" % main_id}
            return json.dumps(info)
        options.update({
            "main_id": str(main_id),
            "task_id": data.task_id,
            "anno_summary_path": anno_summary_path,
            "recode_vcf_path": recode_vcf_path,
            "csv_dir": csv_dir,
            "trait_list": ','.join(trait_list),
            "is_wgs_result": is_wgs_result,
            "update_info": json.dumps(update_info)     # 给work_flow用
        })
        self.set_sheet_data(name="dna_evolution.hapmap", member_id=member_id, project_sn=project_sn,
                            task_id=data.task_id, main_table_name="Hapmap/" + main_table_name, options=options,
                            module_type="module", params=params, db_type="dna_evolution", analysis_name="Hapmap")
        task_info = super(HapmapAction, self).POST()
        task_info['id'] = str(main_id)
        return json.dumps(task_info)

    # def get_target_path(self, client, path):
    #     """
    #     获取远程磁盘的路径
    #     :return:
    #     """
    #     if client not in ['client01', 'client03']:
    #         raise Exception("client{}类型不正确！".format(client))
    #     if client == 'client01':
    #         target_path = os.path.join("/mnt/ilustre/data", path)
    #     else:
    #         target_path = os.path.join("/mnt/ilustre/tsanger-data", path)
    #     return target_path

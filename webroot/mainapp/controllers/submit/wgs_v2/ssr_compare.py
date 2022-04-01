# -*- coding: utf-8 -*-
# __author__: zengjing
# modified: 20190403

import re
import web
import json
import datetime
from types import StringTypes
from collections import defaultdict
from bson.objectid import ObjectId
from mainapp.models.mongo.dna import Dna
from mainapp.libs.signature import check_sig
from mainapp.controllers.project.dna_controller import DnaController


class SsrCompareAction(DnaController):
    """
    SSR比较分析接口
    """
    def __init__(self):
        super(SsrCompareAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        print data
        params = ["task_id", "task_type", "submit_location", "analysis_model", "alle_number", "region",
                  "region_type", "chongmingming_result"]
        for param in params:
            if not hasattr(data, param):
                info = {"success": False, "info": "缺少%s参数" % param}
                return json.dumps(info)
        try:
            allele_num = [int(n) for n in data.alle_number.split(",")]
        except:
            info = {"success": False, "info": "等位基因数目:%s是逗号分隔的正整数" % data.alle_number}
            return json.dumps(info)
        status, info = self.get_region_select(data)
        if not status:
            return info
        new_params = {"allele_num": allele_num, "region_select": info}
        params_json = {
            "task_id": data.task_id,
            "task_type": int(data.task_type),
            "submit_location": data.submit_location,
            "analysis_model": data.analysis_model,
            "alle_number": data.alle_number,
            "region": data.region if data.region_type == "allregion" else json.loads(data.region),
            "region_type": data.region_type,
            "chongmingming_result": data.chongmingming_result
        }
        if hasattr(data, "sample_info"):
            if json.loads(data.sample_info) and json.loads(data.sample_info) != [""]:
                status, info = self.sample_params_check(data)
                if not status:
                    return info
                cs_list, compare_type, genotype, depth = info[0], info[1], info[2], info[3]
                params_json["sample_info"] = json.loads(data.sample_info)
                params_json["sample_genotype"] = json.loads(data.sample_genotype) if data.analysis_model == "multiple" else data.sample_genotype
                params_json["sample_depth"] = json.loads(data.sample_depth) if data.analysis_model == "multiple" else data.sample_depth
                params_json["compare_type"] = json.loads(data.compare_type) if data.analysis_model == "multiple" else data.compare_type
        if data.analysis_model == "single":
            new_params["sample_info"] = cs_list
            new_params["compare_type"] = compare_type
            new_params["gene_type"] = genotype
            new_params["depth"] = depth
        elif data.analysis_model == "multiple":
            if not hasattr(data, "sample_info") and not hasattr(data, "group_dict"):
                info = {"success": False, "info": "多条件组合模式必须选择样本或组进行比较"}
                return json.dumps(info)
            if hasattr(data, "sample_info"):
                if json.loads(data.sample_info) and json.loads(data.sample_info) != [""]:
                    sample_params = defaultdict(dict)
                    for i in range(len(cs_list)):
                        sample_params[cs_list[i]]["compare_type"] = compare_type[i]
                        for j in range(len(cs_list[i].split("_VS_"))):
                            sample_params[cs_list[i]][cs_list[i].split("_VS_")[j]] = {
                                "depth": depth[i].split("|")[j],
                                "gene_type": genotype[i].split("|")[j]
                            }
                    new_params["sample_params"] = sample_params
            if hasattr(data, "group_dict"):
                if json.loads(data.group_dict) and json.loads(data.group_dict) != [""]:
                    status, info = self.group_params_check(data)
                    if not status:
                        return info
                    group_dict, group_depth, group_maf, group_miss = info[0], info[1], info[2], info[3]
                    group_params = {}
                    group_name = data.group_name.split(",")
                    for i in range(len(group_name)):
                        group_params[group_name[i]] = {
                            "sample": group_dict[group_name[i]],
                            "depth": group_depth[i],
                            "gene_fre": group_maf[i],
                            "miss_fre": group_miss[i]
                        }
                    new_params["group_params"] = group_params
                    params_json["group_id"] = json.loads(data.group_id)
                    params_json["group_ids"] = data.group_ids
                    params_json["group_name"] = data.group_name
                    params_json["spe_str"] = data.spe_str
                    params_json["group_dict"] = json.loads(data.group_dict)
                    params_json["group_depth"] = json.loads(data.group_depth)
                    params_json["group_maf"] = json.loads(data.group_maf)
                    params_json["group_miss"] = json.loads(data.group_miss)
        else:
            info = {"success": False, "info": "分析模式analysis_model:%s只能是single或multiple" % data.analysis_model}
            return json.dumps(info)
        params = json.dumps(params_json, sort_keys=True, separators=(',', ':'))
        result = Dna("dna_wgs_v2").find_one(collection="sg_task", query_dic={"task_id": data.task_id})
        try:
            project_sn = result["project_sn"]
            member_id = result["member_id"]
        except:
            info = {"success": False, "info": "sg_task表里没有project_sn, member_id信息，请检查!"}
            return json.dumps(info)
        results = Dna("dna_wgs_v2").find_many(collection="sg_ssr_marker", query_dic={"task_id": data.task_id, "status": "end"})
        for result in results:
            if result["specimen_list"] == ["Reference"]:
                continue
            ssr_vcf = result["ssr_vcf"]
        if data.chongmingming_result:
            main_table_name = data.chongmingming_result
        else:
            main_table_name = "SsrCompare_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        mongo_data = [
            ("name", main_table_name),
            ("status", "start"),
            ("project_sn", project_sn),
            ("task_id", data.task_id),
            ("params", params),
            ("desc", "S比较分析"),
            ("type", data.analysis_model),
            ("created_ts", datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        ]
        main_id = Dna("dna_wgs_v2").insert_main_table(collection="sg_ssr_compare", data=mongo_data)
        Dna("dna_wgs_v2").update_db_record(collection="sg_ssr_compare", query_dict={"_id": main_id},
                                           update_dict={"main_id": main_id})
        update_info = {str(main_id): "sg_ssr_compare"}
        options = {
            "ssr_vcf": ssr_vcf,
            "params_config": json.dumps(new_params, sort_keys=True, separators=(',', ':')),
            "analysis_model": data.analysis_model,
            "update_info": json.dumps(update_info),
            "main_id": str(main_id)
        }
        main_table_name = "SsrCompare/SsrCompare_" + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3])
        self.set_sheet_data(name="wgs_v2.ssr_compare", member_id=member_id, project_sn=project_sn, task_id=data.task_id, main_table_name=main_table_name,
                            options=options, module_type="tool", params=params, db_type="dna_wgs_v2", analysis_name="SsrCompare")
        task_info = super(SsrCompareAction, self).POST()
        task_info['id'] = str(main_id)
        return json.dumps(task_info)

    def get_region_select(self, data):
        """
        对region进行检查，得到参数region_select: {"chr1":["-", "start-end"]}
        """
        if data.region_type not in ["allregion", "location", "custom"]:
            info = {"success": False, "info": "region_type: %s只能是allregion/location/custom" % data.region_type}
            return False, json.dumps(info)
        region_select = defaultdict(list)
        if data.region_type == "allregion":
            region_select = "all"
        elif data.region_type == "location":
            main_id = ObjectId(data.region)
            result = Dna("dna_wgs_v2").find_one(collection="sg_marker_position_table", query_dic={"main_id": main_id})
            if not result:
                info = {"success": False, "info": "sg_marker_position_table表里没有找到id: %s" % main_id}
                return False, json.dumps(info)
            results = Dna("dna_wgs_v2").find(collection="sg_marker_position_table_detail", query_dic={"marker_id": main_id})
            for result in results:
                if not re.search("-", result["location"]):
                    info = {"success": False, "info": "区域:%s 的start和end之间用-分隔" % result["location"]}
                    return False, json.dumps(info)
                chr = result["location"].split(":")[0]
                pos = result["location"].split(":")[1]
                region_select[chr].append(pos)
        else:
            for region in json.loads(data.region):
                chr = region.split(":")[0]
                pos = region.split(":")[1]
                if not re.search("-", pos):
                    info = {"success": False, "info": "区域: %s 的start和end之间用-分隔" % region}
                    return False, json.dumps(info)
                region_select[chr].append(pos)
        return True, region_select

    def sample_params_check(self, data):
        """
        样本间比较参数检查
        """
        single_params = ["sample_info", "sample_genotype", "sample_depth", "compare_type"]
        for param in single_params:
            if not hasattr(data, param):
                info = {"success": False, "info": "缺少样本间比较的%s参数" % param}
                return False, json.dumps(info)
        cs_list = []
        for cs in json.loads(data.sample_info):
            cs_list.append(cs.split("|")[0] + "_VS_" + cs.split("|")[1])
        if data.analysis_model == "multiple":
            compare_type, genotype, depth = [], [], []
            compare_type = json.loads(data.compare_type)
            genotype = json.loads(data.sample_genotype)
            depth = json.loads(data.sample_depth)
            for i in range(len(json.loads(data.sample_info))):
                if compare_type[i] not in ["all", "same", "diff"]:
                    info = {"success": False, "info": "比较类型只能是all/same/diff"}
                    return False, json.dumps(info)
                for g in genotype[i].split("|"):
                    if g not in ["all", "same", "mixed"]:
                        info = {"success": False, "info": "基因型只能是all/same/mixed"}
                        return False, json.dumps(info)
                for dp in depth[i].split("|"):
                    if not re.search("-", dp):
                        info = {"success": False, "info": "测序深度必须用-分隔"}
                        return False, json.dumps(info)
        else:
            compare_type = data.compare_type
            if compare_type not in ["all", "same", "diff"]:
                info = {"success": False, "info": "比较类型只能是all/same/diff"}
                return False, json.dumps(info)
            genotype = data.sample_genotype
            if genotype not in ["all", "same", "mixed"]:
                info = {"success": False, "info": "基因型只能是all/same/mixed"}
                return False, json.dumps(info)
            if not re.search(",", data.sample_depth):
                info = {"success": False, "info": "单条件批量模式的时候测序深度depth:%s必须用逗号分隔" % data.sample_depth}
                return False, info
            dp = data.sample_depth.split(",")
            if dp[0] and dp[1]:
                if int(dp[0]) > int(dp[1]):
                    info = {"success": False, "info": "测序深度depth:%s必须后一个值大于等于前一个值" % dp}
                    return False, info
            depth = data.sample_depth.replace(",", "-")
        return True, [cs_list, compare_type, genotype, depth]

    def group_params_check(self, data):
        """
        组内比较参数检查
        """
        params = ["group_dict", "group_depth", "group_maf", "group_miss", "group_id", "group_ids", "group_name", "spe_str"]
        for param in params:
            if not hasattr(data, param):
                info = {"success": False, "info": "缺少%s参数" % param}
                return False, json.dumps(info)
        group_dict = json.loads(data.group_dict)
        group_depth = json.loads(data.group_depth)
        group_maf = json.loads(data.group_maf)
        group_miss = json.loads(data.group_miss)
        for g in group_dict:
            if not isinstance(group_dict[g], list):
                info = {"success": False, "info": "每个分组都需要是列表" % group_dict[g]}
                return False, json.dumps(info)
        group_name = data.group_name.split(",")
        for i in range(len(group_name)):
            try:
                if not re.search("-", group_depth[i]):
                    info = {"success": False, "info": "深度之间用-分隔" % group_depth[i]}
                    return False, json.dumps(info)
                if not re.search("-", group_maf[i]):
                    info = {"success": False, "info": "基因型频率之间用-分隔" % group_maf[i]}
                    return False, json.dumps(info)
                if float(group_miss[i]) > 1 or float(group_miss[i]) < 0:
                    info = {"success": False, "info": "缺失率:%s 只能是0-1" % group_miss[i]}
                    return False, json.dumps(info)
            except:
                info = {"success": False, "info": "组内比较参数：深度、基因型频率、缺失率列表长度需和分组方案一致"}
                return False, json.dumps(info)
        return True, [group_dict, group_depth, group_maf, group_miss]

# -*- coding: utf-8 -*-
# __author__: zengjing
# modified: 20190225

import re
import web
import json
import datetime
from collections import defaultdict
from bson.objectid import ObjectId
from mainapp.models.mongo.dna import Dna
from mainapp.libs.signature import check_sig
from mainapp.controllers.project.dna_controller import DnaController


class SvCompareAction(DnaController):
    """
    SV比较分析接口
    """
    def __init__(self):
        super(SvCompareAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        print data
        params = ["task_id", "task_type", "submit_location", "analysis_model", "sample_info",
                  "region", "region_type", "genotype", "depth", "compare_type", "chongmingming_result"]
        for param in params:
            if not hasattr(data, param):
                info = {"success": False, "info": "缺少%s参数" % param}
                return json.dumps(info)
        sample_info = json.loads(data.sample_info)
        if len(sample_info) == 0:
            info = {"success": False, "info": "必须选择样本才能进行比较分析"}
            return info
        if data.analysis_model not in ["single", "multiple"]:
            info = {"success": False, "info": "分析模式analysis_model:%s只能是single或multiple" % data.analysis_model}
            return info
        status, region = self.get_region_select(data)
        if not status:
            return region
        if data.analysis_model == "single":
            status, params_option = self.single_option(data)
        else:
            status, params_option = self.multiple_params(data)
        if not status:
            return params_option
        params_option["Region"] = region
        print params_option
        params_json = {
            "task_id": data.task_id,
            "task_type": int(data.task_type),
            "submit_location": data.submit_location,
            "analysis_model": data.analysis_model,
            "sample_info": json.loads(data.sample_info),
            "compare_type": json.loads(data.compare_type) if data.analysis_model == "multiple" else data.compare_type,
            "region_type": data.region_type,
            "genotype": json.loads(data.genotype) if data.analysis_model == "multiple" else data.genotype,
            "depth": json.loads(data.depth) if data.analysis_model == "multiple" else data.depth,
            "chongmingming_result": data.chongmingming_result
        }
        params_json["region"] = "all" if data.region == "all" else json.loads(data.region)
        params = json.dumps(params_json, sort_keys=True, separators=(',', ':'))
        result = Dna("dna_wgs_v2").find_one(collection="sg_task", query_dic={"task_id": data.task_id})
        try:
            project_sn = result["project_sn"]
            member_id = result["member_id"]
            sv_vcf = result["pop_sv_vcf"]
        except:
            info = {"success": False, "info": "sg_task表里没有project_sn, member_id信息，请检查!"}
            return json.dumps(info)
        if data.chongmingming_result:
            main_table_name = data.chongmingming_result
        else:
            main_table_name = "SvCompare_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        mongo_data = [
            ("name", main_table_name),
            ("status", "start"),
            ("project_sn", project_sn),
            ("task_id", data.task_id),
            ("params", params),
            ("type", data.analysis_model),
            ("desc", "S比较分析"),
            ("created_ts", datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        ]
        main_id = Dna("dna_wgs_v2").insert_main_table(collection="sg_sv_compare", data=mongo_data)
        Dna("dna_wgs_v2").update_db_record(collection="sg_sv_compare", query_dict={"_id": main_id},
                                           update_dict={"main_id": main_id})
        update_info = {str(main_id): "sg_sv_compare"}
        options = {
            "sv_vcf": sv_vcf,
            "params_config": json.dumps(params_option),
            "analysis_model": data.analysis_model,
            "update_info": json.dumps(update_info),
            "main_id": str(main_id)
        }
        main_table_name = "SvCompare/SvCompare_" + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3])
        self.set_sheet_data(name="wgs_v2.sv_compare", member_id=member_id, project_sn=project_sn, task_id=data.task_id, main_table_name=main_table_name,
                            options=options, module_type="tool", params=params, db_type="dna_wgs_v2", analysis_name="SvCompare")
        task_info = super(SvCompareAction, self).POST()
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

    def multiple_params(self, data):
        """
        多条件组合模式参数
        """
        sample_info = json.loads(data.sample_info)
        compare_type = json.loads(data.compare_type)
        geno_types = json.loads(data.genotype)
        depths = json.loads(data.depth)
        params_option = {}
        for i in range(len(sample_info)):
            cs_info = sample_info[i].split("|")
            cs_ = cs_info[0] + "_VS_" + cs_info[1]
            if compare_type[i] not in ["same", "diff", "all"]:
                info = {"success": False, "info": "比较类型compare_type:%s只能是same或diff或all" % compare_type[i]}
                return False, info
            params_option[cs_] = {"compare_type": compare_type[i]}
            geno_type = geno_types[i].split("|")
            depth = depths[i].split("|")
            for j in range(len(cs_info)):
                if geno_type[j] not in ["same", "mixed", "all"]:
                    info = {"success": False, "info": "基因型genotype:%s只能是same或mixed或all" % data.genotype}
                    return False, info
                dp = depth[j].split("-")
                if dp[0] and dp[1]:
                    if int(dp[0]) > int(dp[1]):
                        info = {"success": False, "info": "测序深度depth:%s必须后一个值大于等于前一个值" % dp}
                        return False, info
                params_option[cs_][cs_info[j]] = {"geno_type": geno_type[j], "depth": depth[j]}
        return True, params_option

    def single_option(self, data):
        """
        单条件批量模式options检查
        """
        sample_info = json.loads(data.sample_info)
        if data.compare_type not in ["same", "diff", "all"]:
            info = {"success": False, "info": "比较类型compare_type:%s只能是same或diff或all" % data.compare_type}
            return False, info
        if data.genotype not in ["same", "mixed", "all"]:
            info = {"success": False, "info": "基因型genotype:%s只能是same或mixed或all" % data.genotype}
            return False, info
        if not re.search(",", data.depth):
            info = {"success": False, "info": "单条件批量模式的时候测序深度depth:%s必须用逗号分隔" % data.depth}
            return False, info
        dp = data.depth.split(",")
        if dp[0] and dp[1]:
            if int(dp[0]) > int(dp[1]):
                info = {"success": False, "info": "测序深度depth:%s必须后一个值大于等于前一个值" % dp}
                return False, info
        params_option = {}
        for i in range(len(sample_info)):
            cs_info = sample_info[i].split("|")
            cs_ = cs_info[0] + "_VS_" + cs_info[1]
            params_option[cs_] = {"compare_type": data.compare_type}
            for j in range(len(cs_info)):
                params_option[cs_][cs_info[j]] = {"geno_type": data.genotype, "depth": data.depth.replace(",", "-")}
        return True, params_option

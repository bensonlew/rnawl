# -*- coding: utf-8 -*-
# __author__: zengjing
# modified: 20190401

import re
import web
import json
import datetime
from types import StringTypes
from bson.objectid import ObjectId
from collections import defaultdict
from mainapp.models.mongo.dna import Dna
from mainapp.libs.signature import check_sig
from mainapp.controllers.project.dna_controller import DnaController


class RegionAnnoAction(DnaController):
    """
    功能注释接口
    """
    def __init__(self):
        super(RegionAnnoAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        print data
        params = ["task_id", "task_type", "submit_location", "anno_id", "region", "region_type", "chongmingming_result"]
        for param in params:
            if not hasattr(data, param):
                info = {"success": False, "info": "缺少参数：%s,请检查" % param}
                return json.dumps(info)
        if data.anno_id == "all" and data.region_type == "allregion":
            info = {"success": False, "info": "全基因组的注释已经在工作流运行，请查看初始注释结果"}
            return json.dumps(info)
        status, region_select = self.get_region_select(data)
        if not status:
            return region_select
        result = Dna("dna_wgs_v2").find_one(collection="sg_task", query_dic={"task_id": data.task_id})
        try:
            project_sn = result["project_sn"]
            member_id = result["member_id"]
            go_level2 = result["go_summary_path"]
            species_version_id = result["genome_version_id"]
        except:
            info = {"success": False, "info": "sg_task表里没有project_sn, member_id,go_level2,species_version_id信息，请检查!"}
            return json.dumps(info)
        if data.anno_id == "all":
            result = Dna("dna_wgs_v2").find_one(collection="sg_region_anno", query_dic={"task_id": data.task_id, "name": "origin_region_anno"})
        else:
            result = Dna("dna_wgs_v2").find_one(collection="sg_region_anno", query_dic={"main_id": ObjectId(data.anno_id)})
        try:
            pop_summary = result['pop_summary']
        except:
            info = {"success": False, "info": "sg_region_anno表里没有pop_summary: %s信息，请检查!" % data.anno_id}
            return json.dumps(info)
        params_json = {
            "task_id": data.task_id,
            "task_type": int(data.task_type),
            "submit_location": data.submit_location,
            "anno_id": data.anno_id,
            "region_type": data.region_type,
            "chongmingming_result": data.chongmingming_result
        }
        params_json["region"] = "all" if data.region == "all" else json.loads(data.region)
        params = json.dumps(params_json, sort_keys=True, separators=(',', ':'))
        if data.chongmingming_result:
            main_table_name = data.chongmingming_result
        else:
            main_table_name = "RegionAnno_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        mongo_data = [
            ("name", main_table_name),
            ("status", "start"),
            ("project_sn", project_sn),
            ("task_id", data.task_id),
            ("params", params),
            ("desc", "区域功能注释"),
            ("created_ts", datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        ]
        main_id = Dna("dna_wgs_v2").insert_main_table(collection="sg_region_anno", data=mongo_data)
        Dna("dna_wgs_v2").update_db_record(collection="sg_region_anno", query_dict={"_id": main_id},
                                           update_dict={"main_id": main_id})
        update_info = {str(main_id): "sg_region_anno"}
        options = {
            "pop_summary": pop_summary,
            "go_level2": go_level2,
            "region_select": json.dumps(region_select),
            "update_info": json.dumps(update_info),
            "main_id": str(main_id),
            "species_version_id": str(species_version_id)
        }
        main_table_name = "RegionAnno/RegionAnno_" + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3])
        self.set_sheet_data(name="wgs_v2.region_anno", member_id=member_id, project_sn=project_sn, task_id=data.task_id, main_table_name=main_table_name,
                            options=options, module_type="tool", params=params, db_type="dna_wgs_v2", analysis_name="RegionAnno")
        task_info = super(RegionAnnoAction, self).POST()
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

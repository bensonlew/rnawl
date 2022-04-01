# -*- coding: utf-8 -*-
# __author__ = 'qing_mei
# modified 20180816
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


class RegionAnnoAction(DnaController):
    """
    群体进化 基因注释接口
    """

    def __init__(self):
        super(RegionAnnoAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        print data
        params = ["sg_anno_params_id", "chongmingming_result"]
        for param in params:
            if not hasattr(data, param):
                info = {"success": False, "info": "缺少%s参数!" % param}
                return json.dumps(info)
        task_result = Dna("dna_evolution").find_one(collection="sg_task", query_dic={"task_id": data.task_id})
        if not task_result:
            info = {"success": False,
                    "info": "sg_task表里没有task_id: %s对应的信息，请检查!" % data.task_id}
            return json.dumps(info)
        project_sn = task_result["project_sn"]
        member_id = task_result["member_id"]
        pop_summary_path = task_result['pop_summary']   # pop_summary字段先取自于gmap sg_task内，后面确定
        # pop_summary_path = self.get_target_path(data.client, pop_summary_path)
        params_result = Dna("dna_evolution").find_one(collection="sg_anno_params",
                                                      query_dic={"_id": ObjectId(data.sg_anno_params_id)})
        if not params_result:
            info = {"success": False,
                    "info": "sg_anno_params表里没有_id: %s对应的信息，请检查!" % data.sg_anno_params_id}
            return json.dumps(info)
        if "region_path" not in params_result.keys():
            info = {"success": False, "info": "sg_anno_params表里没有region_path，请检查!"}
            return json.dumps(info)
        region_path = params_result["region_path"]
        # region_path = self.get_target_path(data.client, region_path)
        params_json = {
            "sg_anno_params_id": data.sg_anno_params_id,    # 页面另外两个信息前端可自查表，不用存params?!
            "region_select": data.region_select,
            "submit_location": data.submit_location
        }
        params = json.dumps(params_json, sort_keys=True, separators=(',', ':'))     # 返回json
        # main_table_name = 'RegionAnno_' + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")
        main_table_name = Dna("dna_evolution").set_main_table_name("RegionAnno_", data.chongmingming_result)
        mongo_data = [
            ("name", main_table_name),
            ("project_sn", project_sn),
            ("task_id", data.task_id),
            ("created_ts", datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ("status", "start"),
            ("params", params),
            ("desc", "基因注释接口"),
        ]
        main_id = Dna("dna_evolution").insert_main_table(collection="sg_region_anno", data=mongo_data)
        Dna("dna_evolution").update_db_record(collection="sg_region_anno", query_dict={"_id": main_id},
                                              update_dict={"main_id": main_id})
        update_info = {str(main_id): "sg_region_anno"}
        result = Dna("dna_evolution").find_one(collection="sg_region_anno", query_dic={"_id": main_id})
        if not result:
            info = {"success": False,
                    "info": "sg_region_anno表里没有main_id: %s对应的信息，请检查!" % main_id}
            return json.dumps(info)
        task_result = Dna("dna_evolution").find_one(collection="sg_task", query_dic={"task_id": data.task_id})
        if not task_result:
            info = {"success": False,
                    "info": "sg_task表里没有tasl_id: %s对应的信息，请检查!" % data.task_id}
            return json.dumps(info)
        if "genome_version_id" not in task_result.keys():
            info = {"success": False, "info": "sg_task表里没有genome_version_id，请检查!"}
            return json.dumps(info)
        genome_version_id = task_result["genome_version_id"]
        options = {
            "main_id": str(main_id),
            "task_id": data.task_id,
            "pop_summary_path": pop_summary_path,
            "region_path": region_path,
            "region_select": data.region_select,   # ["chr3-1-117536", "chr4-1-118215"]
            "update_info": json.dumps(update_info),     # 给work_flow用
            "genome_version_id": str(genome_version_id)
        }
        self.set_sheet_data(name="dna_evolution.region_anno", member_id=member_id, project_sn=project_sn,
                            task_id=data.task_id, main_table_name="RegionAnno/" + main_table_name, options=options,
                            module_type="module", params=params, db_type="dna_evolution", analysis_name="RegionAnno")
        task_info = super(RegionAnnoAction, self).POST()
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

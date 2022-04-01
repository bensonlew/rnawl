# -*- coding: utf-8 -*-
# __author__: zhaobinbin
# modified: 20200622

import os
import json
import datetime
from api_base import ApiBase
from types import StringTypes
from bson.objectid import ObjectId
from collections import defaultdict


class SsrSpecimen(ApiBase):
    """
    SSR标记导表
    """
    def __init__(self, bind_object):
        super(SsrSpecimen, self).__init__(bind_object)

    def add_sg_ssr_marker(self, project_sn, task_id, params=None, type="Reference"):
        """
        sg_ssr_marker导表
        """
        if not params and type == "Reference":
            params = json.dumps({"specimen_list": ["Reference"], "submit_location": "ssrmarker",
                                 "chongmingming_result": "", "task_id": task_id, "task_type": 2})
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "status": "end",
            "params": params if params else "null",
            "name": type if type else "Reference",
            "desc": "SSR标记",
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "specimen_list":  ["Reference"]
        }
        main_id = self.db["sg_ssr"].insert_one(insert_data).inserted_id
        self.db["sg_ssr"].update({"_id": main_id}, {"$set": {"main_id": main_id}})
        return main_id

    def add_sg_ssr_marker_stat(self, task_id, ssr_id, ssr_stat_path):
        """
        参考基因组的sg_ssr_marker_stat导表
        """
        ssr_id = self.check_objectid(ssr_id)
        self.check_exists(ssr_stat_path)
        type_value = []
        with open(ssr_stat_path, "r") as f:
            lines = f.readlines()
            header = lines[0].strip().split("\t")
            for i in range(1, len(header)):
                type_value.append(0)
            for line in lines[1:]:
                item = line.strip().split("\t")
                for i in range(1, len(header)):
                    type_value[i-1] += int(item[i])
        data_list = []
        name_list = ["c", "c_star", "p1", "p2", "p3", "p4", "p5", "p6"]
        for i in range(8):
            insert_data = {
                "ssr_id": ssr_id,
                "name": name_list[i],
                "value": type_value[i + 1],
                "type": "column",
                "category": ""
            }
            data_list.append(insert_data)
        if len(data_list) == 0:
            self.bind_object.logger.info("{}文件为空！".format(ssr_stat_path))
        else:
            self.col_insert_data("sg_ssr_bar", data_list)
        # type_list = ["c*", "p1", "p2", "p3", "p4", "p5", "p6"] # 可能以后会用得着
        update_dict = {
            "column_data": json.dumps({
                "name": "name",
                "data": "value",
                "category": "category",
                "condition": {'type': "column"}
            })}
        self.update_db_record("sg_ssr", {"_id": ssr_id}, update_dict)
        print "SSR标记导表完成"


if __name__ == "__main__":
    a = SsrSpecimen(None)
    project_sn = "wgs_v2"
    task_id = "sanger_85433"
    ssr_stat_path = "/mnt/ilustre/users/sanger-dev/workspace/20200622/SsrAnalysis_tsg_3421_0622142955635118_1627/output/ref.ssr.stat"
    # ssr_id = a.add_sg_ssr_marker(project_sn, task_id, params=None, type="SsrSpecimen")
    ssr_id = "5cbd7a1e5da1d6a3842a5b12"
    # a.add_sg_ssr_marker_stat(task_id, ssr_id, ssr_stat_path)
    a.add_sg_ssr_marker_stat(task_id, ssr_id, ssr_stat_path)
    # ssr_stat_path = "/mnt/ilustre/users/sanger-dev/workspace/20190402/Single_ssr_compare/SsrCompare/output/BDZ_VS_AA.stat.xls"
    # ssr_detail_path = "/mnt/ilustre/users/sanger-dev/workspace/20190402/Single_ssr_compare/SsrCompare/output/BDZ_VS_AA.detail.xls"
    # name = None
    # compare_id = a.add_sg_ssr_compare(project_sn, task_id, type="multiple")
    # # a.update_info(coll="sg_ssr_compare", query_dict={"main_id": compare_id}, update_dict={"subname": ["BDZ_VS_AA"], "vcf_path": ssr_detail_path})
    # a.update_info(coll="sg_ssr_compare", query_dict={"main_id": compare_id}, update_dict={"vcf_path": ssr_detail_path})
    # a.add_sg_ssr_compare_stat(compare_id, ssr_stat_path, name=name)
    # a.add_sg_ssr_compare_detail(compare_id, ssr_detail_path, name=name)

# -*- coding: utf-8 -*-
# __author__: zengjing
# modified: 20190321

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
        self._project_type = "dna_wgs_v2"
        self.project_sn = self.bind_object.sheet.project_sn
        self.task_id = self.bind_object.sheet.task_id

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
        main_id = self.db["sg_ssr_marker"].insert_one(insert_data).inserted_id
        self.db["sg_ssr_marker"].update({"_id": main_id}, {"$set": {"main_id": main_id}})
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
        data_list = [{
            "ssr_id": ssr_id,
            "sample_id": "Reference",
            "ssr_num": type_value[0],
            "c": type_value[1],
            "c_star": type_value[2],
            "p1": type_value[3],
            "p2": type_value[4],
            "p3": type_value[5],
            "p4": type_value[6],
            "p5": type_value[7],
            "p6": type_value[8]
        }]
        self.db["sg_ssr_marker_stat"].insert_many(data_list)
        specimen_list = ["Reference"]
        type_list = ["c*", "p1", "p2", "p3", "p4", "p5", "p6"]
        bar_id = self.sg_bar(task_id, ssr_id, "", specimen_list, 1, "ssr_type_summation_bar", "", "")
        for i in range(len(type_list)):
            self.sg_bar_detail(bar_id, type_list[i], [type_value[i+1]], "false")
        bar_id = self.sg_bar(task_id, ssr_id, "", header[2:], 1, "ssr_type_histogram_bar", "", "")
        self.sg_bar_detail(bar_id, "Reference", type_value[1:], "false")
        print "SSR标记导表完成"

    def add_sg_ssr_marker_stat_specimen(self, task_id, ssr_id, ssr_stat_path):
        """
        样本的sg_ssr_marker_stat导表
        """
        ssr_id = self.check_objectid(ssr_id)
        self.check_exists(ssr_stat_path)
        specimen_list = []
        specimen_value, type_value = defaultdict(list), defaultdict(list)
        data_list = []
        with open(ssr_stat_path, "r") as f:
            lines = f.readlines()
            header = lines[0].strip().split("\t")
            for line in lines[1:]:
                item = line.strip().split("\t")
                sample_name = item[0]
                specimen_list.append(sample_name)
                specimen_value[sample_name] = [int(i) for i in item[2:]]
                for i in range(2, len(header)):
                    type_value[header[i]].append(int(item[i]))
                insert_data = {
                    "ssr_id": ssr_id,
                    "sample_id": item[0],
                    "ssr_num": int(item[1]),
                    "c": int(item[2]),
                    "c_star": int(item[3]),
                    "p1": int(item[4]),
                    "p2": int(item[5]),
                    "p3": int(item[6]),
                    "p4": int(item[7]),
                    "p5": int(item[8]),
                    "p6": int(item[9])
                }
                data_list.append(insert_data)
        if data_list:
            self.db["sg_ssr_marker_stat"].insert_many(data_list)
            bar_id = self.sg_bar(task_id, ssr_id, "", specimen_list, 1, "ssr_type_summation_bar", "", "")
            for t in type_value:
                self.sg_bar_detail(bar_id, t, type_value[t], "false")
            bar_id = self.sg_bar(task_id, ssr_id, "", header[2:], 1, "ssr_type_histogram_bar", "", "")
            for s in specimen_value:
                self.sg_bar_detail(bar_id, s, specimen_value[s], "false")
            self.bind_object.logger.info("SSR标记导表完成")
        else:
            self.bind_object.logger.info("SSR标记结果为空")

    def add_sg_ssr_compare(self, project_sn, task_id, type, params=None):
        """
        sg_ssr_compare导表
        """
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "status": "end",
            "params": params if params else "null",
            "type": type,
            "name": "SsrCompare",
            "desc": "SSR比较分析",
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        }
        main_id = self.db["sg_ssr_compare"].insert_one(insert_data).inserted_id
        self.db["sg_ssr_compare"].update({"_id": main_id}, {"$set": {"main_id": main_id}})
        return main_id

    def update_info(self, coll, query_dict, update_dict):
        """
        更新sg_ssr_compare表的single下的subname、vcf_path，或multiple下的header、header_dict、vcf_path
        """
        self.db[coll].update(query_dict, {"$set": update_dict})

    def add_sg_ssr_compare_stat(self, compare_id, ssr_stat_path, name=None):
        """
        sg_ssr_compare_stat导表
        有name的时候是单条件批量模式，没有的时候是多条件组合模式
        """
        compare_id = self.check_objectid(compare_id)
        self.check_exists(ssr_stat_path)
        data_list = []
        with open(ssr_stat_path, "r") as f:
            lines = f.readlines()
            header = lines[0].strip().split("\t")
            for line in lines[1:]:
                item = line.strip().split("\t")
                insert_data = {
                    "compare_id": compare_id,
                    "chr": item[0],
                    "ssr_num": int(item[1]),
                    "c": int(item[2]),
                    "c_star": int(item[3]),
                    "p1": int(item[4]),
                    "p2": int(item[5]),
                    "p3": int(item[6]),
                    "p4": int(item[7]),
                    "p5": int(item[8]),
                    "p6": int(item[9])
                }
                if name:
                    insert_data["name"] = name
                data_list.append(insert_data)
        if data_list:
            self.db["sg_ssr_compare_stat"].insert_many(data_list)
        print "SSR比较分析统计导表完成"

    def add_sg_ssr_compare_detail(self, compare_id, ssr_detail_path, name=None):
        """
        sg_ssr_compare_detail导表
        有name的时候是单条件批量模式，没有的时候是多条件组合模式
        """
        compare_id = self.check_objectid(compare_id)
        self.check_exists(ssr_detail_path)
        data_list, header_list, header_dict = [], [], {}
        with open(ssr_detail_path, "r") as f:
            lines = f.readlines()
            header = lines[0].strip().split("\t")
            for line in lines[1:]:
                item = line.strip().split("\t")
                insert_data = {
                    "compare_id": compare_id,
                    "name": name,
                    "chr": item[0],
                    "ssr_id": item[3],
                    "start": int(item[1]),
                    "end": int(item[2]),
                    "region": item[0] + ":" + item[1] + "-" + item[2],
                    "repeat_unit": item[4],
                    "ref_count": int(item[5])
                }
                for i in range(6, len(header)):
                    insert_data[header[i]] = item[i]
                if name:
                    insert_data["name"] = name
                data_list.append(insert_data)
        if data_list:
            self.db["sg_ssr_compare_detail"].insert_many(data_list)
        else:
            self.bind_object.logger.info("SSR比较结果为空")
        if not name:
            for i in header[6:]:
                if i.endswith("_repeat_count"):
                    header_list.append(i.split("_repeat_count")[0] + " Repeat count")
                    header_dict[i] = i.split("_repeat_count")[0] + " Repeat count"
                elif i.endswith("_allele_depth"):
                    header_list.append(i.split("_allele_depth")[0] + " Allele Depth")
                    header_dict[i] = i.split("_allele_depth")[0] + " Allele Depth"
                elif i.endswith("_allele_frequency"):
                    header_list.append(i.split("_allele_frequency")[0] + " Allele Frequency")
                    header_dict[i] = i.split("_allele_frequency")[0] + " Allele Frequency"
            self.db["sg_ssr_compare"].update({"main_id": compare_id}, {"$set": {"header": header_list, "header_dict": header_dict}})
            self.bind_object.logger.info("SSR比较header更新完成")


if __name__ == "__main__":
    a = SsrSpecimen(None)
    project_sn = "wgs_v2"
    task_id = "sanger_85433"
    ssr_stat_path = "/mnt/ilustre/users/sanger-dev/workspace/20190409/WgsV2_sanger_85433/output/08.ssr/ssr.stat"
    # ssr_id = a.add_sg_ssr_marker(project_sn, task_id, params=None, type="SsrSpecimen")
    ssr_id = "5cbd7a1e5da1d6a3842a5b12"
    # a.add_sg_ssr_marker_stat(task_id, ssr_id, ssr_stat_path)
    a.add_sg_ssr_marker_stat_specimen(task_id, ssr_id, ssr_stat_path)
    # ssr_stat_path = "/mnt/ilustre/users/sanger-dev/workspace/20190402/Single_ssr_compare/SsrCompare/output/BDZ_VS_AA.stat.xls"
    # ssr_detail_path = "/mnt/ilustre/users/sanger-dev/workspace/20190402/Single_ssr_compare/SsrCompare/output/BDZ_VS_AA.detail.xls"
    # name = None
    # compare_id = a.add_sg_ssr_compare(project_sn, task_id, type="multiple")
    # # a.update_info(coll="sg_ssr_compare", query_dict={"main_id": compare_id}, update_dict={"subname": ["BDZ_VS_AA"], "vcf_path": ssr_detail_path})
    # a.update_info(coll="sg_ssr_compare", query_dict={"main_id": compare_id}, update_dict={"vcf_path": ssr_detail_path})
    # a.add_sg_ssr_compare_stat(compare_id, ssr_stat_path, name=name)
    # a.add_sg_ssr_compare_detail(compare_id, ssr_detail_path, name=name)

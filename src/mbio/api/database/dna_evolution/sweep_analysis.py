# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# last_modify:20180917

from api_base import ApiBase
from collections import defaultdict
import os
import datetime
import json


class SweepAnalysis(ApiBase):
    """
    选择性消除导表
    """
    def __init__(self, bind_object):
        super(SweepAnalysis, self).__init__(bind_object)
        self._project_type = "dna_evolution"
        self.project_sn = self.bind_object.sheet.project_sn
        self.task_id = self.bind_object.sheet.id.split("_SweepAnalysis")[0]

    def add_sg_sweep(self, diff_group, params=None, name=None):
        """
        主表sg_sweep
        差异分组列表，如["Q1_vs_Q3"]
        """
        if params:
            params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        else:
            params = None
        insert_data = {
            "project_sn": self.project_sn,
            "task_id": self.task_id,
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "status": "end",
            "params": params if params else "",
            "name": name if name else "origin_sweep_analysis",
            "diff_group": diff_group,
            "desc": "选择性消除",
        }
        main_id = self.db['sg_sweep'].insert_one(insert_data).inserted_id
        self.update_db_record("sg_sweep", {"_id": main_id}, {"main_id": main_id})
        return main_id

    def add_sg_sweep_detail(self, sweep_id, pi_tajimad_fst_path, variant_num_path, diff_group_name, window_step):
        """
        sg_sweep_detail
        pi_tajimad_fst_path:1-2.pi_tajimaD_fst.detail
        variant_num_path:vcf.variant.txt
        diff_group_name:差异分组，A_vs_B
        window_step: 滑窗步长
        """
        sweep_id = self.check_objectid(sweep_id)
        self.check_exists(pi_tajimad_fst_path)
        self.check_exists(variant_num_path)
        varinat_info, manhan_data = {}, {}
        data_list = []
        with open(variant_num_path, "r") as f:
            lines = f.readlines()
            for line in lines[1:]:
                item = line.strip().split("\t")
                chr = item[0].replace('"', "")
                varinat_info[chr + "_" + item[1]] = int(item[2])
        pi1 = diff_group_name.split("_vs_")[0] + "_pi"
        pi2 = diff_group_name.split("_vs_")[1] + "_pi"
        tajima1 = diff_group_name.split("_vs_")[0] + "_tajima"
        tajima2 = diff_group_name.split("_vs_")[1] + "_tajima"
        axis_type = [pi1, pi2, tajima1, tajima2, "fst", "θπ"]
        combina_data = {pi1: [], pi2: [], tajima1: [], tajima2: [], "fst": [], "θπ": []}
        switch = [tajima1 + "/" + pi1, tajima2 + "/" + pi2, "θπ/fst"]
        chr_list = []
        with open(pi_tajimad_fst_path, "r") as f:
            lines = f.readlines()
            for line in lines[1:]:
                item = line.strip().split(" ")
                chr = item[0].replace('"', "")
                fst = round(float(item[-1]), 4)
                pi1_data = round(float(item[2]), 4)
                pi2_data = round(float(item[3]), 4)
                tajima1_data = round(float(item[4]), 4)
                tajima2_data = round(float(item[5]), 4)
                the = round(float(item[2])/float(item[3]), 4) if float(item[3]) != 0 else 0
                insert_data = {
                    "sweep_id": sweep_id,
                    "diff_group_name": diff_group_name,
                    "chr": chr,
                    "start": int(item[1]),
                    "end": int(item[1]) + window_step,
                    "fst": fst,
                    pi1: pi1_data,
                    pi2: pi2_data,
                    tajima1: tajima1_data,
                    tajima2: tajima2_data
                }
                try:
                    insert_data["variant_num"] = varinat_info[chr + "_" + item[1]]
                except:
                    insert_data["variant_num"] = 0
                data_list.append(insert_data)
                if chr not in chr_list:
                    chr_list.append(chr)
                    manhan_data[chr] = {"pos_data": [], "pi1": [], "pi2": [], "tajima1": [], "tajima2": [], "fst": []}
                manhan_data[chr]["pos_data"].append(int(item[1]))
                manhan_data[chr]["pi1"].append(pi1_data)
                manhan_data[chr]["pi2"].append(pi2_data)
                manhan_data[chr]["tajima1"].append(tajima1_data)
                manhan_data[chr]["tajima2"].append(tajima2_data)
                manhan_data[chr]["fst"].append(fst)
                combina_data[pi1].append(pi1_data)
                combina_data[pi2].append(pi2_data)
                combina_data[tajima1].append(tajima1_data)
                combina_data[tajima2].append(tajima2_data)
                combina_data["fst"].append(fst)
                combina_data["θπ"].append(the)
        if data_list:
            self.col_insert_data("sg_sweep_detail", data_list)
            combination_id = self.add_sg_combination(sweep_id, diff_group_name, switch)
            self.add_sg_combination_detail(combination_id, combina_data, axis_type)
            # self.add_sg_manhattan_detail(sweep_id, diff_group_name, chr_list, manhan_data)
        else:
            self.bind_object.logger.info("{}文件结果为空".format(pi_tajimad_fst_path))

    def add_sg_manhattan(self, sweep_id, name, ext_name, chr_list):
        """
        sg_manhattan
        name:一个比较组
        """
        sweep_id = self.check_objectid(sweep_id)
        insert_data = {
            "origin_id": sweep_id,
            "task_id": self.task_id,
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "name": name,
            "type": 1,
            "location" : "sweep_manhattan",
            "ext_name": ext_name,
            "chr_list": chr_list
        }
        main_id = self.db['sg_manhattan'].insert_one(insert_data).inserted_id
        return main_id

    def add_sg_manhattan_detail(self, sweep_id, name, chr_list, manhan_data):
        """
        sg_manhattan_detail
        name: pi/tajima/fst
        """
        sweep_id = self.check_objectid(sweep_id)
        data_list = []
        for ext_name in ["pi1", "pi2", "tajima1", "tajima2", "fst"]:
            manhattan_id = self.add_sg_manhattan(sweep_id, name, ext_name, chr_list)
            for chr in chr_list:
                insert_data = {
                    "manhattan_id": manhattan_id,
                    "name": chr,
                    "pos_data": manhan_data[chr]["pos_data"],
                    "value": manhan_data[chr][ext_name]
                }
                data_list.append(insert_data)
        self.col_insert_data("sg_manhattan_detail", data_list)

    def add_sg_manhattan_path(self, sweep_id, name, ext_name, png_path, pdf_path):
        """
        sg_manhattan,修改成静态图
        name:一个比较组
        ext_name:比较组中的pop1或pop2
        """
        sweep_id = self.check_objectid(sweep_id)
        insert_data = {
            "origin_id": sweep_id,
            "task_id": self.task_id,
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "name": name,
            "type": 1,
            "location" : "sweep_manhattan",
            "ext_name": ext_name,
            "png_path": png_path,
            "pdf_path": pdf_path
        }
        main_id = self.db['sg_manhattan'].insert_one(insert_data).inserted_id

    def add_sg_combination(self, sweep_id, diff_group_name, axis_type):
        """
        sg_combination
        diff_group_name:一个比较组
        """
        sweep_id = self.check_objectid(sweep_id)
        insert_data = {
            "origin_id": sweep_id,
            "task_id": self.task_id,
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "name": diff_group_name,
            "type": 1,
            "axis_type": axis_type,
            "location" : "sweep"
        }
        main_id = self.db['sg_combination'].insert_one(insert_data).inserted_id
        return main_id

    def add_sg_combination_detail(self, combination_id, combina_data, axis_type):
        """
        sg_combination_detail
        """
        combination_id = self.check_objectid(combination_id)
        data_list = []
        for axis in axis_type:
            insert_data = {
                "combination_id": combination_id,
                "name": axis,
                "value": combina_data[axis]
            }
            data_list.append(insert_data)
        self.col_insert_data("sg_combination_detail", data_list)

    def update_sweep_path(self, sweep_id, sweep_dir, vcf_path):
        """
        更新sweep_dir,vcf_path路径
        """
        sweep_id = self.check_objectid(sweep_id)
        # self.check_exists(sweep_dir)  # 对象存储的时候，这个方法获取不到文件路径
        # self.check_exists(vcf_path)
        self.update_db_record("sg_sweep", {"main_id": sweep_id}, {"sweep_dir": sweep_dir, "vcf_path": vcf_path})


if __name__ == "__main__":
    a = SweepAnalysis(None)
    project_sn = "test_zj"
    task_id = "tsg_32120"
    diff_group_name = "1_vs_2"
    window_step = 10000
    pi_tajimad_fst_path = "/mnt/ilustre/users/sanger-dev/workspace/20181009/Single_sweep_analysis2/SweepAnalysis/output/1-2.pi_tajimaD_fst.detail"
    variant_num_path = "/mnt/ilustre/users/sanger-dev/workspace/20181009/Single_sweep_analysis2/SweepAnalysis/vcf.variant.txt"
    sweep_id = a.add_sg_sweep(project_sn, task_id)
    a.add_sg_sweep_detail(task_id, sweep_id, pi_tajimad_fst_path, variant_num_path, diff_group_name, window_step)

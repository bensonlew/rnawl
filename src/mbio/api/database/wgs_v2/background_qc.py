# -*- coding: utf-8 -*-
# __author__: zengjing
# modified: 20190305

import os
import json
import datetime
from types import StringTypes
from bson.objectid import ObjectId
from biocluster.api.database.base import Base, report_check
from mbio.api.database.wgs.api_base import ApiBase


class BackgroundQc(ApiBase):
    """
    WGS V2项目背景和质控导表
    """
    def __init__(self, bind_object):
        super(BackgroundQc, self).__init__(bind_object)
        self._project_type = "dna_wgs_v2"
        self.project_sn = ''
        self.task_id = ''

    def add_sg_results(self, project_sn, task_id, path):
        """
        sg_results导表
        """
        data_list = []
        with open(path, "r") as f:
            lines = f.readlines()
            for line in lines[1:]:
                item = line.strip().split("\t")
                insert_data = {
                    "project_sn": project_sn,
                    "task_id": task_id,
                    "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                    "total_data": item[0],
                    "average_data": item[1],
                    "average_q30": item[2],
                    "average_depth": item[3],
                    "mapping_efficiency": item[4],
                    "average_coverage": item[5],
                    "average_gene": item[6],
                    "total_snp": int(item[7]),
                    "total_indel": int(item[8]),
                    "average_snp": int(item[9]),
                    "average_indel": int(item[10]),
                    "average_sv": int(item[11]),
                    "average_cnv": int(item[12])
                }
                data_list.append(insert_data)
        self.col_insert_data(collection="sg_results", data_list=data_list)

    def add_sg_specimen(self):
        """
        sg_specimen导表
        """
        self.project_sn = self.bind_object.sheet.project_sn
        self.task_id = self.bind_object.sheet.id
        result = self.col_find_one("sg_task", {"task_id": self.task_id})
        if not result:
            raise Exception("没有在sg_task表中找到task_id:{}".format(self.task_id))
        results = self.col_find("sg_specimen_other", {"task_id": self.task_id, "selected": "1"})
        if not results:
            raise Exception("没有在sg_specimen_other表中找到task_id:{}".format(self.task_id))
        data_list = []
        for result in results:
            insert_data = {
                "project_sn": self.project_sn,
                "task_id": self.task_id,
                "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                "analysis_name": result["analysis_name"],
                "initial_name": result["initial_name"],
                "new_name": result["analysis_name"],
                "library": result["library"],
                "batch": result["batch"],
                "raw_data": result["file_name"],
                "desc": ""
            }
            data_list.append(insert_data)
        self.col_insert_data(collection="sg_specimen", data_list=data_list)

    def add_sg_software(self, project_sn, task_id, path):
        """
        sg_software导表
        """
        data_list = []
        with open(path, "r") as f:
            lines = f.readlines()
            for line in lines[1:]:
                item = line.strip().split("\t")
                insert_data = {
                    "project_sn": project_sn,
                    "task_id": task_id,
                    "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                    "item": str(item[0]),
                    "software_name": item[1],
                    "version": item[2],
                    "params": item[3]
                }
                print insert_data
                data_list.append(insert_data)
        self.col_insert_data(collection="sg_software", data_list=data_list)

    def add_fastp_json_stat(self, project_sn, task_id, json_dir):
        """
        fastp进行质控后直接用json文件进行导表
        """
        for f in os.listdir(json_dir):
            specimen_id = f.split(".json")[0]
            sample_id = specimen_id.split("-")[0]
            json_path = os.path.join(json_dir, f)
            r = open(json_path, "r")
            json_dict = json.loads(r.read())
            summary = json_dict["summary"]
            raw_stat = summary["before_filtering"]
            clean_stat = summary["after_filtering"]
            query_dict = {"task_id": task_id, "specimen_id": sample_id}
            result = self.col_find_one("sg_specimen_qc", query_dict)
            raw_gc_content = round(float(raw_stat["gc_content"]) * 100, 2)
            raw_q30_rate = round(float(raw_stat["q30_rate"]) * 100, 2)
            clean_gc_rate = round(float(clean_stat["gc_content"]) * 100, 2)
            clean_q30_rate = round(float(clean_stat["q30_rate"]) * 100, 2)
            if result:
                main_id = result["_id"]
                raw_reads = int(raw_stat["total_reads"]) + result["raw_reads"]
                raw_base = int(raw_stat["total_bases"]) + result["raw_base"]
                raw_gc_rate = round((raw_gc_content + result["raw_gc_rate"]) / 2, 2)
                raw_q30_rate = round((raw_q30_rate + result["raw_q30_rate"]) / 2, 2)
                clean_reads = int(clean_stat["total_reads"]) + result["clean_reads"]
                clean_base = int(clean_stat["total_bases"]) + result["clean_base"]
                clean_gc_rate = round((clean_gc_rate + result["clean_gc_rate"]) / 2, 2)
                clean_q30_rate = round((clean_q30_rate + result["clean_q30_rate"]) / 2, 2)
                update_dict = {"raw_reads": raw_reads, "raw_base": raw_base, "raw_gc_rate": raw_gc_rate,
                               "raw_q30_rate": raw_q30_rate,
                               "clean_reads": clean_reads, "clean_base": clean_base, "clean_gc_rate": clean_gc_rate,
                               "clean_q30_rate": clean_q30_rate, "main_id": main_id}
                self.update_db_record(collection="sg_specimen_qc", query_dict=query_dict, update_dict=update_dict)
            else:
                insert_data = {
                    "project_sn": project_sn,
                    "task_id": task_id,
                    "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                    "specimen_id": sample_id,
                    "raw_reads": int(raw_stat["total_reads"]),
                    "raw_base": int(raw_stat["total_bases"]),
                    "raw_gc_rate": raw_gc_content,
                    "raw_q30_rate": raw_q30_rate,
                    "clean_reads": int(clean_stat["total_reads"]),
                    "clean_base": int(clean_stat["total_bases"]),
                    "clean_gc_rate": clean_gc_rate,
                    "clean_q30_rate": clean_q30_rate
                }
                main_id = self.db["sg_specimen_qc"].insert_one(insert_data).inserted_id
                self.db["sg_specimen_qc"].update_one(query_dict, {"$set": {"main_id": main_id}})
            raw_read1 = json_dict["read1_before_filtering"]
            raw_read1_cate = raw_read1["total_cycles"]
            raw_read1_content = raw_read1["content_curves"]
            raw_read1_a = raw_read1_content["A"]
            raw_read1_t = raw_read1_content["T"]
            raw_read1_c = raw_read1_content["C"]
            raw_read1_g = raw_read1_content["G"]
            raw_read1_n = raw_read1_content["N"]
            raw_read1_mean = raw_read1["quality_curves"]["mean"]
            raw_read2 = json_dict["read2_before_filtering"]
            raw_read2_cate = raw_read2["total_cycles"]
            raw_read2_content = raw_read2["content_curves"]
            categories, raw_e_list = [], []
            for i in range(1, raw_read1_cate):
                categories.append(i)
            for i in range(raw_read1_cate, raw_read1_cate + raw_read2_cate):
                categories.append(i)
            categories.append(raw_read1_cate + raw_read2_cate)
            raw_read1_a.extend(raw_read2_content["A"])
            raw_read1_t.extend(raw_read2_content["T"])
            raw_read1_c.extend(raw_read2_content["C"])
            raw_read1_g.extend(raw_read2_content["G"])
            raw_read1_n.extend(raw_read2_content["N"])
            raw_read1_a = self.return_new_list(raw_read1_a)
            raw_read1_t = self.return_new_list(raw_read1_t)
            raw_read1_c = self.return_new_list(raw_read1_c)
            raw_read1_g = self.return_new_list(raw_read1_g)
            raw_read1_n = self.return_new_list(raw_read1_n)
            raw_read1_mean.extend(raw_read2["quality_curves"]["mean"])
            for mean in raw_read1_mean:
                raw_e_list.append(10 ** (float(mean)/(-10)) * 100)
            try:
                ext_title = "-" + specimen_id.split("-")[1]
            except:
                ext_title = ""
            curve_id = self.sg_curve(task_id, main_id, sample_id, categories, 1, "raw_reads", '',
                                     {"type": "sample", "title": sample_id, "ext_title": ext_title})
            self.update_db_record("sg_curve", {"_id": curve_id}, {"x_data": raw_read1_cate + raw_read2_cate})
            self.sg_curve_detail(curve_id, "A", raw_read1_a)
            self.sg_curve_detail(curve_id, "T", raw_read1_t)
            self.sg_curve_detail(curve_id, "C", raw_read1_c)
            self.sg_curve_detail(curve_id, "G", raw_read1_g)
            self.sg_curve_detail(curve_id, "N", raw_read1_n)
            curve_id = self.sg_curve(task_id, main_id, sample_id, categories, 1, "error_raw_reads", '',
                                     {"type": "sample", "title": sample_id, "ext_title": ext_title})
            self.update_db_record("sg_curve", {"_id": curve_id}, {"x_data": raw_read1_cate + raw_read2_cate})
            self.sg_curve_detail(curve_id, "error", raw_e_list)
            clean_read1 = json_dict["read1_after_filtering"]
            clean_read1_cate = clean_read1["total_cycles"]
            clean_read1_content = clean_read1["content_curves"]
            clean_read1_a = clean_read1_content["A"]
            clean_read1_t = clean_read1_content["T"]
            clean_read1_c = clean_read1_content["C"]
            clean_read1_g = clean_read1_content["G"]
            clean_read1_n = clean_read1_content["N"]
            clean_read2 = json_dict["read2_after_filtering"]
            clean_read2_cate = clean_read2["total_cycles"]
            clean_read2_content = clean_read2["content_curves"]
            clean_read1_mean = clean_read1["quality_curves"]["mean"]
            clean_read1_a.extend(clean_read2_content["A"])
            clean_read1_t.extend(clean_read2_content["T"])
            clean_read1_c.extend(clean_read2_content["C"])
            clean_read1_g.extend(clean_read2_content["G"])
            clean_read1_n.extend(clean_read2_content["N"])
            clean_read1_a = self.return_new_list(clean_read1_a)
            clean_read1_t = self.return_new_list(clean_read1_t)
            clean_read1_c = self.return_new_list(clean_read1_c)
            clean_read1_g = self.return_new_list(clean_read1_g)
            clean_read1_n = self.return_new_list(clean_read1_n)
            clean_read1_mean.extend(clean_read2["quality_curves"]["mean"])
            categories, clean_e_list = [], []
            for mean in clean_read1_mean:
                clean_e_list.append(10 ** (float(mean)/(-10)) * 100)
            for i in range(1, clean_read1_cate):
                categories.append(i)
            for i in range(clean_read1_cate, clean_read1_cate + clean_read2_cate):
                categories.append(i)
            categories.append(clean_read1_cate + clean_read2_cate)
            curve_id = self.sg_curve(task_id, main_id, sample_id, categories, 1, "clean_reads", '',
                                     {"type": "sample", "title": sample_id, "ext_title": ext_title})
            self.update_db_record("sg_curve", {"_id": curve_id}, {"x_data": raw_read1_cate + raw_read2_cate})
            self.sg_curve_detail(curve_id, "A", clean_read1_a)
            self.sg_curve_detail(curve_id, "T", clean_read1_t)
            self.sg_curve_detail(curve_id, "C", clean_read1_c)
            self.sg_curve_detail(curve_id, "G", clean_read1_g)
            self.sg_curve_detail(curve_id, "N", clean_read1_n)
            curve_id = self.sg_curve(task_id, main_id, sample_id, categories, 1, "error_clean_reads", '',
                                     {"type": "sample", "title": sample_id, "ext_title": ext_title})
            self.update_db_record("sg_curve", {"_id": curve_id}, {"x_data": raw_read1_cate + raw_read2_cate})
            self.sg_curve_detail(curve_id, "error", clean_e_list)

    def return_new_list(self, old_list):
        """
        将列表里的值乘以100再返回
        """
        return [i * 100 for i in old_list]


if __name__ == "__main__":
    a = BackgroundQc(None)
    # a.add_sg_results(project_sn="wgs_v2", task_id="wgs_v2-1", path="/mnt/ilustre/users/sanger-dev/sg-users/xuanhongdong/wgs/api_data/sg_results.txt")
    # a.add_sg_specimen(project_sn="wgs_v2", task_id="wgs_v2-1", path="/mnt/ilustre/users/sanger-dev/sg-users/xuanhongdong/wgs/api_data/cleandata/fqlist.txt")
    # a.add_sg_software(project_sn="wgs_v2", task_id="wgs_v2-1", path="/mnt/ilustre/users/sanger-dev/sg-users/xuanhongdong/wgs/api_data/sg_software.txt")
    # a.add_fastp_json_stat(project_sn="wgs_v2", task_id="wgs_v2-1", json_dir="/mnt/ilustre/users/sanger-dev/workspace/20190305/Single_qc_stat/QcStat/qc_stat")
    project_sn = "wgs_v2"
    task_id = "sanger_85433"
    json_dir = "/mnt/ilustre/users/sanger-dev/workspace/20190409/WgsV2_sanger_85433/QcStat/output/qc_stat"
    a.add_fastp_json_stat(project_sn, task_id, json_dir)

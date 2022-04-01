# -*- coding: utf-8 -*-
# __author__ = "zengjing"
# last_modify: 20180624

import os
import re
import json
import types
import datetime
from bson.objectid import ObjectId
from collections import defaultdict
from biocluster.config import Config
from mbio.files.medical.html import HtmlFile
from biocluster.api.database.base import Base, report_check

class DatasplitNew(Base):
    def __init__(self, bind_object):
        super(DatasplitNew, self).__init__(bind_object)
        self._project_type = "datasplit"

    def check_objectid(self, id_):
        """
        用于检查并转成成ObjectID
        """
        if not isinstance(id_, ObjectId):
            if isinstance(id_, types.StringTypes):
                id_ = ObjectId(id_)
            else:
                raise Exception("id必须为ObjectId对象或其对应的字符串!")
        return id_

    def check_exists(self, path):
        """
        检查文件是否存在
        """
        if not os.path.exists(path):
            raise Exception("{}所指定的路径不存在，请检查".format(path))
        else:
            return True

    # @report_check
    def add_flowcell_summary(self, split_id, lane_html, lane_barcode_html):
        """
        导入文库统计信息
        lane_html:lane.html
        lane_barcode_html:laneBarcode.html
        """
        self.check_exists(lane_html)
        self.check_exists(lane_barcode_html)
        parser = HtmlFile()
        parser.set_path(lane_html)
        parser.get_info()
        tab_list = parser.tab_list
        split_id = self.check_objectid(split_id)
        insert_data = {
            "split_id": split_id,
            "clusters_raw": tab_list[1][1][0],
            "clusters_pf": tab_list[1][1][1],
            "yield_nbases": tab_list[1][1][2]
        }
        self.db["sg_split_summary"].insert_one(insert_data)
        self.bind_object.logger.info("导入flowcell_summary成功")
        data_list = []
        for info in tab_list[2][1:]:
            insert_data = {
                "split_id": split_id,
                "lane": info[0],
                "clusters_pf": info[1],
                "lane_rate": info[2],
                "perfect_barcode_rate": info[3],
                "mis_barcode": info[4],
                "yield": info[5],
                "clusters_pf_rate": info[6],
                "base_q30": info[7],
                "quality_score": info[8]
            }
            data_list.append(insert_data)
        self.db["sg_split_lane_summary"].insert_many(data_list)
        self.bind_object.logger.info("导入lane_summary成功")
        parser2 = HtmlFile()
        parser2.set_path(lane_barcode_html)
        parser2.get_info()
        tab_list2 = parser2.tab_list
        data_list = []
        for info in tab_list2[2][1:]:
            insert_data = {
                "split_id": split_id,
                "lane": info[0],
                "project": info[1],
                "library_name": info[2],
                "barcode_seq": info[3],
                "clusters_pf": info[4],
                "lane_rate": info[5]
            }
            if len(info) == 7:
                insert_data["perfect_barcode"] = ""
                insert_data["mis_barcode"] = ""
                insert_data["yield"] = info[6]
                insert_data["clusters_pf_rate"] = ""
                insert_data["base_q30"] = ""
                insert_data["quality_score"] = ""
            else:
                insert_data["perfect_barcode"] = info[6]
                insert_data["mis_barcode"] = info[7]
                insert_data["yield"] = info[8]
                insert_data["clusters_pf_rate"] = info[9]
                insert_data["base_q30"] = info[10]
                insert_data["quality_score"] = info[11]
            data_list.append(insert_data)
        self.db["sg_split_lane_summary_detail"].insert_many(data_list)
        self.bind_object.logger.info("导入lane_summary_detail成功")
        data_list = []  # 第四张表格有点特殊，第一列合并lane,需要进行一些处理
        head = tab_list2[3][1]
        lane_list = []
        for i in range(len(head) / 3):
            lane_list.append(head[i*3])
        for info in tab_list2[3][1:]:
            if len(info) != 0:
                if len(info) == len(head):
                    for i in range(len(lane_list)):
                        lane = lane_list[i]
                        insert_data = {
                            "split_id": split_id,
                            "lane": lane,
                            "count": info[i*3+1],
                            "squence": info[i*3+2]
                        }
                        data_list.append(insert_data)
                else:
                    for i in range(len(lane_list)):
                        lane = lane_list[i]
                        insert_data = {
                            "split_id": split_id,
                            "lane": lane,
                            "count": info[i*2],
                            "squence": info[i*2+1]
                        }
                        data_list.append(insert_data)
        self.db["sg_split_unknow_barcode"].insert_many(data_list)
        self.bind_object.logger.info("导入top_unknown_barcodesl成功")
        self.bind_object.logger.info("一次拆分结果导入成功")

    def update_lib_path(self, split_id, fastq_dir, sample_sheet, s3_upload_dir):
        """
        一拆拆分后更新拆分出来的文库路径
        fastq_dir: 一拆拆分后所有文库存放的文件夹路径
        sample_sheet: 进行拆分的拆分表
        s3_upload_dir:上传到s3对象存储的dir
        """
        split_id = self.check_objectid(split_id)
        self.check_exists(fastq_dir)
        self.check_exists(sample_sheet)
        sample_sheet_dict = {}
        with open(sample_sheet, "rb") as f:
            lines = f.readlines()
            for line in lines[1:]:
                item = line.strip().split(",")
                sample_sheet_dict[item[1]] = item[2]
        for lib in os.listdir(fastq_dir):
            lib_id = sample_sheet_dict[lib]
            lib_dir = os.path.join(fastq_dir, lib)
            md5sum_file = os.path.join(lib_dir,"md5sum.txt")
            md5sum_dict = {}
            if os.path.exists(md5sum_file):
                with open(md5sum_file,"r") as mf:
                    while 1:
                        line = mf.readline()
                        if not line:
                            break
                        fd = line.rstrip().split("  ")
                        md5sum_dict[fd[1]] = fd[0]
            fq_path, work_path,md5info = [], [],[]
            for fq in os.listdir(lib_dir):
                if s3_upload_dir:
                    path = os.path.join(s3_upload_dir, lib + "/" + fq)
                else:
                    path = os.path.join(lib_dir, lib + "/" + fq)
                path_ = os.path.join(lib_dir, fq)
                if re.search(r"R1.raw", fq):
                    md5info.insert(0,md5sum_dict[fq])
                    fq_path.insert(0, path)
                    work_path.insert(0, path_)
                elif re.search(r"_R1_", fq):
                    md5info.insert(0,md5sum_dict[fq])
                    fq_path.insert(0, path)
                    work_path.insert(0, path_)
                elif re.search(r"md5sum", fq):
                    continue
                else:
                    md5info.append(md5sum_dict[fq])
                    fq_path.append(path)
                    work_path.append(path_)
            query_dict = {"split_id": split_id, "library_number": lib_id}
            # print query_dict
            # print fq_path
            self.db["sg_split_library"].update(query_dict, {"$set": {"path": ";".join(fq_path), "work_path": ";".join(work_path),"md5sum":";".join(md5info)}})
            query_dict = {"split_id": split_id, "library_number": lib_id}
            spe_result = self.db["sg_split_specimen"].find(query_dict)
            if spe_result.count() > 0:
                raw_bytes = []
                for f_ in work_path:
                    raw_bytes.append(str(os.path.getsize(f_)))
                # if spe_result.count() == 1 and spe_result[0]["product_type"] not in ["meta", "dna"]:
                if spe_result.count() == 1:
                    self.db["sg_split_specimen"].update(query_dict, {"$set": {"raw_md5sum":";".join(md5info),"raw_path": ";".join(fq_path), "work_path": ";".join(work_path), "raw_bytes": ";".join(raw_bytes)}})
                # if spe_result[0]["product_type"] == "meta":
                #     self.db["sg_split_specimen"].update(query_dict, {"$set": {"raw_path": ";".join(fq_path), "work_path": ";".join(work_path), "raw_bytes": ";".join(raw_bytes)}},False,True)
        self.bind_object.logger.info("文库路径更新成功!")

    def update_sample_path(self, split_id, fastq_dir, s3_upload_dir, product_type):
        """
        更新样本的路径
        """
        split_id = self.check_objectid(split_id)
        self.check_exists(fastq_dir)
        md5sum_file = os.path.join(fastq_dir,"md5sum.txt")
        md5sum_dict = {}
        if os.path.exists(md5sum_file):
            with open(md5sum_file,"r") as m5:
                while 1:
                    line = m5.readline()
                    if not line:
                        break 
                    fd = line.rstrip().split("  ")
                    md5sum_dict[fd[1]] = fd[0]
        if product_type == "meta_raw":
            sample_info, work_path,md5sum_info = defaultdict(list), defaultdict(list),defaultdict(list)
            for f in os.listdir(fastq_dir):
                if not f.endswith(".gz"):
                    continue
                fq_path = s3_upload_dir+f if s3_upload_dir.endswith("/") else s3_upload_dir+"/"+f
                sample_id = ObjectId(f.split("--")[2])
                if f.endswith(".R1.raw.fastq.gz"):
                    md5sum_info[sample_id].insert(0, md5sum_dict[f])
                    sample_info[sample_id].insert(0, fq_path)
                    work_path[sample_id].insert(0, os.path.join(fastq_dir, f))
                elif f.endswith(".R2.raw.fastq.gz"):
                    md5sum_info[sample_id].append(md5sum_dict[f])
                    sample_info[sample_id].append(fq_path)
                    work_path[sample_id].append(os.path.join(fastq_dir, f))
            for sample_id in sample_info:
                raw_bytes = []
                for f_ in work_path[sample_id]:
                    raw_bytes.append(str(os.path.getsize(f_)))
                update_dict = {"raw_md5sum":";".join(md5sum_info[sample_id]),"raw_path": ";".join(sample_info[sample_id]), "work_path": ";".join(work_path[sample_id]), "raw_bytes": ";".join(raw_bytes)}
                self.db["sg_split_specimen"].update({"_id": sample_id}, {"$set": update_dict})
        elif product_type == "meta_clean":
            for f in os.listdir(fastq_dir):
                if f.endswith("fastq.gz"):
                    sample_id = ObjectId(f.split("--")[2])
                    s3_path = s3_upload_dir+f if s3_upload_dir.endswith("/") else s3_upload_dir+"/"+f
                    clean_bytes = str(os.path.getsize(os.path.join(fastq_dir, f)))
                    self.db["sg_split_specimen"].update({"_id": sample_id}, {"$set": {"clean_md5sum":md5sum_dict[f],"clean_path": s3_path, "clean_bytes": clean_bytes, "clean_work_path": os.path.join(fastq_dir, f)}})
        elif product_type == "dna_raw":
            sample_info, work_path,md5sum_info = defaultdict(list), defaultdict(list),defaultdict(list)
            for f in os.listdir(fastq_dir):
                fq_path = s3_upload_dir+f if s3_upload_dir.endswith("/") else s3_upload_dir+"/"+f
                if f.endswith(".R1.raw.fastq.gz"):
                    lib_spe = f.split(".R1.raw.fastq.gz")[0]
                    md5sum_info[lib_spe].insert(0,md5sum_dict[f])
                    sample_info[lib_spe].insert(0, fq_path)
                    work_path[lib_spe].insert(0, os.path.join(fastq_dir, f))
                elif f.endswith(".R2.raw.fastq.gz"):
                    lib_spe = f.split(".R2.raw.fastq.gz")[0]
                    md5sum_info[lib_spe].append(md5sum_dict[f])
                    sample_info[lib_spe].append(fq_path)
                    work_path[lib_spe].append(os.path.join(fastq_dir, f))
            for lib_spe in sample_info:
                raw_bytes = []
                for f_ in work_path[lib_spe]:
                    raw_bytes.append(str(os.path.getsize(f_)))
                update_dict = {"raw_md5sum":";".join(md5sum_info[lib_spe]),"raw_path": ";".join(sample_info[lib_spe]), "work_path": ";".join(work_path[lib_spe]), "raw_bytes": ";".join(raw_bytes)}
                if len(lib_spe.split("--")) > 2:
                    sample_id = lib_spe.split("--")[0]
                    lib = lib_spe.split("--")[1]
                    sample = lib_spe.split("--")[2]
                    self.db["sg_split_specimen"].update({"_id": ObjectId(sample_id)}, {"$set": update_dict})
                else:
                    lib = lib_spe.split("--")[0]
                    sample = lib_spe.split("--")[1]
                    self.db["sg_split_specimen"].update({"split_id": split_id, "library_number": lib, "specimen_name": sample}, {"$set": update_dict})
        elif product_type in ["microbial_genome", "meta_genomic", "rna", "prokaryotic_rna", "lncrna", "dna"]:
            sample_info, sample_work_info, sample_size,md5sum_info = defaultdict(list), defaultdict(list), defaultdict(list), defaultdict(list)
            for f in os.listdir(fastq_dir):
                fq_path = s3_upload_dir+f if s3_upload_dir.endswith("/") else s3_upload_dir+"/"+f
                fq_work_path = os.path.join(fastq_dir, f)
                if f.endswith(".clean.1.fastq.gz"):
                    lib_spe = f.split(".clean.1.fastq.gz")[0]
                    md5sum_info[lib_spe].insert(0, md5sum_dict[f])
                    sample_info[lib_spe].insert(0, fq_path)
                    sample_work_info[lib_spe].insert(0, fq_work_path)
                    sample_size[lib_spe].insert(0, str(os.path.getsize(os.path.join(fastq_dir, f))))
                elif f.endswith(".clean.2.fastq.gz"):
                    lib_spe = f.split(".clean.2.fastq.gz")[0]
                    md5sum_info[lib_spe].append(md5sum_dict[f])
                    sample_info[lib_spe].append(fq_path)
                    sample_work_info[lib_spe].append(fq_work_path)
                    sample_size[lib_spe].append(str(os.path.getsize(os.path.join(fastq_dir, f))))
            for lib_spe in sample_info:
                lib = lib_spe.split("--")[0]
                sample = lib_spe.split("--")[1]
                query_dict = {"split_id": split_id, "library_number": lib,"specimen_name": sample}
                if product_type == "dna" and len(lib_spe.split("--")) == 3:
                    query_dict = {"split_id": split_id, "library_number": lib, "specimen_name": sample, "order_sn": lib_spe.split("--")[2]}
                update_dict = {"clean_md5sum":";".join(md5sum_info[lib_spe]),"clean_path": ";".join(sample_info[lib_spe]), "clean_work_path": ";".join(sample_work_info[lib_spe]), "clean_bytes": ";".join(sample_size[lib_spe])}
                self.db["sg_split_specimen"].update(query_dict, {"$set": update_dict})
        elif product_type == "mirna":
            for f in os.listdir(fastq_dir):
                if f.endswith(".fasta.gz"):
                    lib = f.split("--")[0]
                    sample = f.split("--")[1].split(".fasta.gz")[0]
                    fa_path = s3_upload_dir+f if s3_upload_dir.endswith("/") else s3_upload_dir+"/"+f
                    fa_work_path = os.path.join(fastq_dir, f)
                    clean_bytes = str(os.path.getsize(os.path.join(fastq_dir, f)))
                    query_dict = {"split_id": split_id, "library_number": lib, "specimen_name": sample}
                    update_dict = {"clean_md5sum":md5sum_dict[f],"clean_path": fa_path, "clean_work_path": fa_work_path, "clean_bytes": clean_bytes}
                    self.db["sg_split_specimen"].update(query_dict, {"$set": update_dict})
                if f.endswith(".fq.gz"):
                    lib = f.split("--")[0]
                    sample = f.split("--")[1].split(".fq.gz")[0]
                    raw75_bytes = str(os.path.getsize(os.path.join(fastq_dir, f)))
                    fa_path = s3_upload_dir+f if s3_upload_dir.endswith("/") else s3_upload_dir+"/"+f
                    query_dict = {"split_id": split_id, "library_number": lib, "specimen_name": sample}
                    update_dict = {"raw75_path": fa_path, "raw75_bytes": raw75_bytes,"raw75_md5sum": md5sum_dict[f]}
                    self.db["sg_split_specimen"].update(query_dict, {"$set": update_dict})

    def update_raw_sample_path(self, split_id):
        """
        根据split_id,将sg_split_library的path更新到sg_split_specimen的raw_path
        解决从二次拆分开始，但二次拆分直接跳过，开始进行样本质控的时候sg_split_specimen的raw_path为空的问题
        """
        split_id = self.check_objectid(split_id)
        result = self.db["sg_split"].find_one({"_id": split_id})
        if result["split_type"] == "first_split" or result["split_type"] == "second_split":
            results = self.db["sg_split_library"].find({"split_id": split_id})
            for result in results:
                results1 = self.db["sg_split_specimen"].find({"split_id": split_id, "library_id": result["_id"]})
                # for result1 in results1:
                #     if result1["raw_path"] == "" and result1["product_type"] != "meta":
                #         self.db["sg_split_specimen"].update({"_id": result1["_id"]}, {"$set": {"raw_path": result["path"], "work_path": result["work_path"]}})
                if results1.count() == 1:
                    result1 = results1[0]
                    if result1["raw_path"] == "":
                        self.db["sg_split_specimen"].update({"_id": result1["_id"]}, {"$set": {"raw_path": result["path"], "work_path": result["work_path"]}})
        self.bind_object.logger.info("更新原始样本路径成功")

    def add_sg_split_clean_qc(self, split_id, fastq_stat, raw_fastq_stat):
        """
        多样性质控后样本统计导表
        """
        split_id = self.check_objectid(split_id)
        self.check_exists(fastq_stat)
        self.check_exists(raw_fastq_stat)
        data_list, sample_ids, raw_data_list = [], [], []
        with open(fastq_stat, "rb") as f:
            lines = f.readlines()
            for line in lines[1:]:
                item = line.strip().split("\t")
                sample = item[0].split("--")
                if sample[2] in sample_ids:
                    continue
                sample_ids.append(sample[2])
                spe_result = self.db["sg_split_specimen"].find_one({"_id": ObjectId(sample[2])})
                # clean_data = round(float(item[2]) / 1014 / 1024, 4)
                # clean_data = round(float(item[2]) / 1000 / 1000, 4)
                clean_data = self.friendly_size(float(item[2]))
                insert_data = {
                    "split_id": split_id,
                    "library_number": sample[1],
                    "specimen_name": sample[3],
                    "project_sn": sample[0],
                    "product_type": "meta",
                    "insert_len": spe_result["insert_size"],
                    "specimen_id": spe_result["_id"],
                    "seq_model": "PE",
                    "total_reads": int(item[1]),
                    "total_bases": int(item[2]),
                    # "clean_data": str(clean_data) + "M",
                    "clean_data": clean_data,
                    "a_rate": float(item[5]),
                    "t_rate": float(item[6]),
                    "c_rate": float(item[7]),
                    "g_rate": float(item[8]),
                    "n_rate": float(item[9]),
                    "gc_rate": float(item[13]),
                    "q20_rate": float(item[11]),
                    "q30_rate": float(item[12]),
                    "error_rate": float(item[10])
                }
                data_list.append(insert_data)
        with open(raw_fastq_stat, "rb") as f:
            lines = f.readlines()
            for line in lines[1:]:
                item = line.strip().split("\t")
                sample = item[0].split("--")
                # if sample[2] in sample_ids:
                #     continue
                # sample_ids.append(sample[2])
                spe_result = self.db["sg_split_specimen"].find_one({"_id": ObjectId(sample[2])})
                clean_data = self.friendly_size(float(item[2]))
                insert_data = {
                    "split_id": split_id,
                    "library_number": sample[1],
                    "specimen_name": sample[3],
                    "project_sn": sample[0],
                    "product_type": "meta",
                    "insert_len": spe_result["insert_size"],
                    "specimen_id": spe_result["_id"],
                    "seq_model": "PE",
                    "total_reads": int(item[1]),
                    "total_bases": int(item[2]),
                    "raw_data": clean_data,
                    "a_rate": float(item[5]),
                    "t_rate": float(item[6]),
                    "c_rate": float(item[7]),
                    "g_rate": float(item[8]),
                    "n_rate": float(item[9]),
                    "gc_rate": float(item[13]),
                    "q20_rate": float(item[11]),
                    "q30_rate": float(item[12]),
                    "error_rate": float(item[10])
                }
                raw_data_list.append(insert_data)
        spe_results = self.db["sg_split_specimen"].find({"split_id": split_id, "product_type" : "meta"})
        for spe_result in spe_results:
            sample_id = str(spe_result["_id"])
            if sample_id not in sample_ids:
                insert_data = {
                    "split_id": split_id,
                    "library_number": spe_result["library_number"],
                    "specimen_name": spe_result["specimen_name"],
                    "project_sn": spe_result["project_sn"],
                    "product_type": "meta",
                    "insert_len": spe_result["insert_size"],
                    "specimen_id": spe_result["_id"],
                    "seq_model": "PE",
                    "total_reads": 0,
                    "total_bases": 0,
                    "clean_data": 0,
                    "a_rate": 0,
                    "t_rate": 0,
                    "c_rate": 0,
                    "g_rate": 0,
                    "n_rate": 0,
                    "gc_rate": 0,
                    "q20_rate": 0,
                    "q30_rate": 0,
                    "error_rate": 0
                }
                data_list.append(insert_data)
        self.db["sg_split_clean_qc"].insert_many(data_list)
        if len(raw_data_list) > 0:
            self.db["sg_split_raw_qc"].insert_many(raw_data_list)
        self.bind_object.logger.info("多样性质控统计导表成功")

    def add_sg_split_qc_mirna(self, split_id, fastq_stat, sample_qc_stat):
        """
        mirna质控样本统计导表
        """
        split_id = self.check_objectid(split_id)
        self.check_exists(fastq_stat)
        self.check_exists(sample_qc_stat)
        raw_data_list, clean_data_list = [], []
        with open(fastq_stat, "rb") as f:
            lines = f.readlines()
            for line in lines[1:]:
                item = line.strip().split("\t")
                lib_num = item[0].split("--")[0]
                sample_name = item[0].split("--")[1]
                spe_result = self.db["sg_split_specimen"].find_one({"split_id": split_id, "library_number": lib_num, "specimen_name": sample_name})
                # raw_data = round(float(item[2]) / 1024 / 1024, 4)
                # raw_data = round(float(item[2]) / 1000 / 1000, 4)
                raw_data = self.friendly_size(float(item[2]))
                insert_data = {
                    "split_id": split_id,
                    "library_number": lib_num,
                    "specimen_name": sample_name,
                    "project_sn": spe_result["project_sn"],
                    "product_type": "mirna",
                    "insert_len": spe_result["insert_size"],
                    "seq_model": "SE",
                    "total_reads": int(item[1]),
                    "total_bases": int(item[2]),
                    # "raw_data": str(raw_data) + "M",
                    "raw_data": raw_data,
                    "a_rate": float(item[5]),
                    "t_rate": float(item[6]),
                    "c_rate": float(item[7]),
                    "g_rate": float(item[8]),
                    "n_rate": float(item[9]),
                    "gc_rate": float(item[13]),
                    "q20_rate": float(item[11]),
                    "q30_rate": float(item[12]),
                    "error_rate": float(item[10]),
                }
                raw_data_list.append(insert_data)
            self.db["sg_split_raw_qc"].insert_many(raw_data_list)
        with open(sample_qc_stat, "rb") as f:
            lines = f.readlines()
            for line in lines[1:]:
                item = line.strip().split("\t")
                lib_num = item[0].split("--")[0]
                sample_name = item[0].split("--")[1]
                spe_result = self.db["sg_split_specimen"].find_one({"split_id": split_id, "library_number": lib_num, "specimen_name": sample_name})
                # raw_data = round(float(item[2]) / 1024 / 1024, 4)
                # raw_data = round(float(item[2]) / 1000 / 1000, 4)
                raw_data = self.friendly_size(float(item[2]))
                insert_data = {
                    "split_id": split_id,
                    "library_number": lib_num,
                    "specimen_name": sample_name,
                    "project_sn": spe_result["project_sn"],
                    "product_type": "mirna",
                    "insert_len": spe_result["insert_size"],
                    "seq_model": "SE",
                    "specimen_id": spe_result["_id"],
                    "total_reads": int(item[6]),
                    "adapter": float(item[7]),
                    "adapter_only": int(item[2]),
                    "18nt": int(item[4]),
                    "32nt": int(item[5])
                }
                clean_data_list.append(insert_data)
            self.db["sg_split_clean_qc"].insert_many(clean_data_list)
        self.bind_object.logger.info("mirna质控统计导表成功")

    def add_mirna_sg_bar(self, split_id, fasta_length_dir):
        """
        mirna长度分布图导表
        """
        split_id = self.check_objectid(split_id)
        self.check_exists(fasta_length_dir)
        for f in os.listdir(fasta_length_dir):
            fasta_length = os.path.join(fasta_length_dir, f)
            lib_num = f.split(".fasta_length.xls")[0].split("--")[0]
            sample_name = f.split(".fasta_length.xls")[0].split("--")[1]
            spe_result = self.db["sg_split_clean_qc"].find_one({"split_id": split_id, "library_number": lib_num, "specimen_name": sample_name})
            name = sample_name
            origin_id = spe_result["_id"]
            location = "mirnabar"
            categories, value_list = [], []
            with open(fasta_length, "rb") as r:
                lines = r.readlines()
                for line in lines[1:]:
                    item = line.strip().split("\t")
                    categories.append(item[0])
                    value_list.append(int(item[1]))
            bar_id = self.add_sg_bar(split_id, origin_id, name, location, categories, types=1)
            self.add_sg_bar_detail(bar_id, name, value_list)
        self.bind_object.logger.info("mirna柱状图导表成功")

    def add_fastp_json_stat(self, split_id, json_dir, product_type):
        """
        fastp进行质控后直接用json文件进行导表
        """
        split_id = self.check_objectid(split_id)
        self.check_exists(json_dir)
        for f in os.listdir(json_dir):
            s_info = f.split(".json")[0].split("--")
            lib = s_info[0]
            sample_name = s_info[1]
            query_dict = {"split_id": split_id, "library_number": lib, "specimen_name": sample_name}
            if product_type == "dna" and len(s_info) == 3:
                query_dict = {"split_id": split_id, "library_number": lib, "specimen_name": sample_name, "order_sn": s_info[2]}
            json_path = os.path.join(json_dir, f)
            r = open(json_path, "r")
            json_dict = json.loads(r.read())
            stat_num = self.fastp_stat(json_dict, err_rate=0.05)
            summary = json_dict["summary"]
            raw_stat = summary["before_filtering"]
            clean_stat = summary["after_filtering"]
            adapter_reads = int(json_dict["adapter_cutting"]["adapter_trimmed_reads"])/2
            # clean_data = round(float(raw_stat["total_bases"]) / 1024 / 1024, 4)
            # clean_data = round(float(raw_stat["total_bases"]) / 1000 / 1000, 4)
            raw_data = self.friendly_size(float(raw_stat["total_bases"]))
            spe_result = self.db["sg_split_specimen"].find_one(query_dict)
            raw_insert = {
                "split_id": split_id,
                "specimen_id": spe_result["_id"],
                "library_number": spe_result["library_number"],
                "specimen_name": spe_result["specimen_name"],
                "project_sn": spe_result["project_sn"],
                "product_type": product_type,
                "insert_len": spe_result["insert_size"],
                "seq_model": "PE",
                # "raw_data": str(clean_data) + "M",
                "raw_data": raw_data,
                "total_reads": int(raw_stat["total_reads"]) / 2,
                "total_bases": int(raw_stat["total_bases"]),
                "gc_rate": float(raw_stat["gc_content"]) * 100,
                "q20_rate": float(raw_stat["q20_rate"]) * 100,
                "q30_rate": float(raw_stat["q30_rate"]) * 100,
                "at_seprate": stat_num["raw_at"],
                "gc_septate": stat_num["raw_gc"],
                "err_num": stat_num["raw_q"],
            }
            # clean_data = round(float(clean_stat["total_bases"]) / 1014 / 1024 / 1024, 4)
            # clean_data = round(float(clean_stat["total_bases"]) / 1000 / 1000 / 1000, 4)
            clean_data = self.friendly_size(float(clean_stat["total_bases"]))
            clean_insert = {
                "split_id": split_id,
                "specimen_id": spe_result["_id"],
                "library_number": spe_result["library_number"],
                "specimen_name": spe_result["specimen_name"],
                "project_sn": spe_result["project_sn"],
                "product_type": product_type,
                "insert_len": spe_result["insert_size"],
                "seq_model": "PE",
                "total_reads": int(clean_stat["total_reads"]) / 2,
                "total_bases": int(clean_stat["total_bases"]),
                # "clean_data": str(clean_data) + "G",
                "clean_data": clean_data,
                "gc_rate": float(clean_stat["gc_content"]) * 100,
                "q20_rate": float(clean_stat["q20_rate"]) * 100,
                "q30_rate": float(clean_stat["q30_rate"]) * 100,
                "n_rate": round(float(json_dict["filtering_result"]["too_many_N_reads"])/int(clean_stat["total_reads"]), 4) * 100,
                "dup_rate": float(json_dict["duplication"]["rate"]) * 100,
                "at_seprate": stat_num["clean_at"],
                "gc_septate": stat_num["clean_gc"],
                "err_num": stat_num["clean_q"],
                "adapter_only": adapter_reads,
                "adapter": round(adapter_reads/(float(raw_stat["total_reads"])/2), 4) * 100
            }
            raw_id = self.db["sg_split_raw_qc"].insert_one(raw_insert).inserted_id
            clean_id = self.db["sg_split_clean_qc"].insert_one(clean_insert).inserted_id
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
            raw_read1_mean.extend(raw_read2["quality_curves"]["mean"])
            for mean in raw_read1_mean:
                raw_e_list.append(10 ** (float(mean)/(-10)))
            curve_id = self.add_sg_curve(raw_id, sample_name, "raw_reads", categories, 1)
            self.add_sg_curve_detail(curve_id, "A", raw_read1_a)
            self.add_sg_curve_detail(curve_id, "T", raw_read1_t)
            self.add_sg_curve_detail(curve_id, "C", raw_read1_c)
            self.add_sg_curve_detail(curve_id, "G", raw_read1_g)
            self.add_sg_curve_detail(curve_id, "N", raw_read1_n)
            curve_id = self.add_sg_curve(raw_id, sample_name, "error_raw_reads", categories, 1)
            self.add_sg_curve_detail(curve_id, sample_name, raw_e_list)
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
            clean_read1_mean.extend(clean_read2["quality_curves"]["mean"])
            categories, clean_e_list = [], []
            for mean in clean_read1_mean:
                clean_e_list.append(10 ** (float(mean)/(-10)))
            for i in range(1, clean_read1_cate):
                categories.append(i)
            for i in range(clean_read1_cate, clean_read1_cate + clean_read2_cate):
                categories.append(i)
            categories.append(clean_read1_cate + clean_read2_cate)
            curve_id = self.add_sg_curve(clean_id, sample_name, "clean_reads", categories, 1)
            self.add_sg_curve_detail(curve_id, "A", clean_read1_a)
            self.add_sg_curve_detail(curve_id, "T", clean_read1_t)
            self.add_sg_curve_detail(curve_id, "C", clean_read1_c)
            self.add_sg_curve_detail(curve_id, "G", clean_read1_g)
            self.add_sg_curve_detail(curve_id, "N", clean_read1_n)
            curve_id = self.add_sg_curve(clean_id, sample_name, "error_clean_reads", categories, 1)
            self.add_sg_curve_detail(curve_id, sample_name, clean_e_list)

    def fastp_stat(self, json_dict, err_rate=0.05):
        """
        统计fastp的A-T/G-C分离率和碱基错误率
        err_rate: A-T/G-C
        """
        raw_read1 = json_dict["read1_before_filtering"]["content_curves"]
        raw_read1_a, raw_read1_t = raw_read1["A"], raw_read1["T"]
        raw_read1_c, raw_read1_g = raw_read1["C"], raw_read1["G"]
        raw_read2 = json_dict["read2_before_filtering"]["content_curves"]
        raw_read2_a, raw_read2_t = raw_read2["A"], raw_read2["T"]
        raw_read2_c, raw_read2_g = raw_read2["C"], raw_read2["G"]
        raw_read1_q = json_dict["read1_before_filtering"]["quality_curves"]["mean"]
        raw_read2_q = json_dict["read2_before_filtering"]["quality_curves"]["mean"]
        clean_read1 = json_dict["read1_after_filtering"]["content_curves"]
        clean_read1_a, clean_read1_t = clean_read1["A"], clean_read1["T"]
        clean_read1_c, clean_read1_g = clean_read1["C"], clean_read1["G"]
        clean_read2 = json_dict["read2_after_filtering"]["content_curves"]
        clean_read2_a, clean_read2_t = clean_read2["A"], clean_read2["T"]
        clean_read2_c, clean_read2_g = clean_read2["C"], clean_read2["G"]
        clean_read1_q = json_dict["read1_before_filtering"]["quality_curves"]["mean"]
        clean_read2_q = json_dict["read2_before_filtering"]["quality_curves"]["mean"]
        raw1_at_num, raw1_gc_num = self.get_list_err_num(raw_read1_a, raw_read1_t, raw_read1_c, raw_read1_g, err_rate=err_rate, start=10)
        raw2_at_num, raw2_gc_num = self.get_list_err_num(raw_read2_a, raw_read2_t, raw_read2_c, raw_read2_g, err_rate=err_rate, start=10)
        clean1_at_num, clean1_gc_num = self.get_list_err_num(clean_read1_a, clean_read1_t, clean_read1_c, clean_read1_g, err_rate=err_rate, start=10)
        clean2_at_num, clean2_gc_num = self.get_list_err_num(clean_read2_a, clean_read2_t, clean_read2_c, clean_read2_g, err_rate=err_rate, start=10)
        raw1_q_num = self.get_list_q20_num(raw_read1_q, q_limit=20, start=10)
        raw2_q_num = self.get_list_q20_num(raw_read2_q, q_limit=20, start=10)
        clean1_q_num = self.get_list_q20_num(clean_read1_q, q_limit=20, start=10)
        clean2_q_num = self.get_list_q20_num(clean_read2_q, q_limit=20, start=10)
        stat_num = {
            "raw_at": raw1_at_num + raw2_at_num,
            "raw_gc": raw1_gc_num + raw2_gc_num,
            "raw_q": raw1_q_num + raw2_q_num,
            "clean_at": clean1_at_num + clean2_at_num,
            "clean_gc": clean1_gc_num + clean2_gc_num,
            "clean_q": clean1_q_num + clean2_q_num
        }
        return stat_num

    def add_meta_sg_bar(self, split_id, trim_hist_dir):
        """
        meta文库长度分布图导表
        """
        split_id = self.check_objectid(split_id)
        self.check_exists(trim_hist_dir)
        lib_files = {}
        for f in os.listdir(trim_hist_dir):
            trim_hist = os.path.join(trim_hist_dir, f)
            lib_num = f.split("--")[1]
            if lib_num not in lib_files.keys():
                lib_files[lib_num] = []
            lib_files[lib_num].append(trim_hist)
        for lib_num in lib_files.keys():
            lib_result = self.db["sg_split_library"].find_one({"split_id": split_id, "library_number": lib_num})
            name = lib_num
            origin_id = lib_result["_id"]
            location = "metabar"
            categories_list, value_dict = [], {}
            for trim_hist in lib_files[lib_num]:
                with open(trim_hist, "rb") as r:
                    lines = r.readlines()
                    for line in lines:
                        item = line.strip().split("\t")
                        categories_list.append(int(item[0]))
                        value_dict[int(item[0])] = int(item[1])
            categories_list = list(set(categories_list))
            categories_list.sort()
            categories, value_list = [], []
            for i in categories_list:
                categories.append(str(i))
                value_list.append(value_dict[i])
            bar_id = self.add_sg_bar(split_id, origin_id, name, location, categories, types=1)
            self.add_sg_bar_detail(bar_id, name, value_list)
        self.bind_object.logger.info("meta长度分布柱状图导表成功")

    def get_list_q20_num(self, q_list, q_limit=20, start=10):
        """
        求列表中大于err_rate的值的数目
        """
        num = 0
        for i in range(start, len(q_list)):
            if float(q_list[i]) < 20:
                num += 1
        return num

    def get_list_err_num(self, a_list, t_list, c_list, g_list, err_rate, start=10):
        """
        求列表中大于err_rate的值的数目
        """
        at_num, gc_num = 0, 0
        for i in range(start, len(a_list)):
            at = abs(float(a_list[i]) - float(t_list[i]))
            gc = abs(float(g_list[i]) - float(c_list[i]))
            if at > err_rate:
                at_num += 1
            if gc > err_rate:
                gc_num += 1
        return at_num, gc_num

    def add_sg_curve(self, origin_id, name, location, categories, types=1):
        """
        导入曲线图数据
        origin_id: sg_split_raw_qc/sg_split_clean_qc的_id
        name: 样本名称
        """
        if not isinstance(origin_id, ObjectId):
            if isinstance(origin_id, StringTypes):
                origin_id = ObjectId(origin_id)
            else:
                raise Exception("origin_id必须为ObjectId对象或其对应的字符串!")
        insert_data = {
            "task_id": "",
            "origin_id": origin_id,
            "name": name,
            "categories": categories,
            "type": types,
            "location": location,
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        }
        return self.db["sg_curve"].insert_one(insert_data).inserted_id

    def add_sg_curve_detail(self, curve_id, name, value):
        """
        导入曲线图细节表
        """
        data_list = []
        new_value = [i * 100 for i in value]  # fastp的结果乘以100得到百分比
        insert_data = {
            "curve_id": curve_id,
            "name": name,
            "value": new_value
        }
        data_list.append(insert_data)
        self.db["sg_curve_detail"].insert_many(data_list)

    def add_sg_bar(self, split_id, origin_id, name, location, categories, types=1):
        """
        导入柱状图数据
        origin_id: sg_split_clean_qc的_id
        name: 样本名称
        """
        if not isinstance(origin_id, ObjectId):
            if isinstance(origin_id, StringTypes):
                origin_id = ObjectId(origin_id)
            else:
                raise Exception("origin_id必须为ObjectId对象或其对应的字符串!")
        insert_data = {
            "task_id": split_id,
            "origin_id": origin_id,
            "name": name,
            "categories": categories,
            "type": types,
            "location": location,
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        }
        return self.db["sg_bar"].insert_one(insert_data).inserted_id

    def add_sg_bar_detail(self, bar_id, name, value):
        """
        导入柱状图细节表
        """
        data_list = []
        insert_data = {
            "bar_id": bar_id,
            "name": name,
            "value": value
        }
        data_list.append(insert_data)
        self.db["sg_bar_detail"].insert_many(data_list)

    def update_cpc_info(self, split_id, output_dir):
        """
        更新sg_split_clean_qc表的比对nt和rfam信息
        """
        split_id = self.check_objectid(split_id)
        self.check_exists(output_dir)
        for f in os.listdir(output_dir):
            f_ = os.path.join(output_dir, f)
            update_dict = {}
            if re.search(r"rfam", f):
                specimen_id = self.check_objectid(f.split("--rfam_summary.xls")[0])
                with open(f_, "rb") as r:
                    for line in r:
                        if re.search(r"rRNA", line):
                            item = line.strip().split("\t")
                            update_dict = {"rRNA": item[1] + "(" + item[2] + ")"}
            elif re.search(r"nt", f):
                specimen_id = self.check_objectid(f.split("--nt_species_stat.xls")[0])
                species_list, num_list, percent_list = [], [], []
                with open(f_, "rb") as r:
                    lines = r.readlines()
                    for line in lines[1:]:
                        item = line.strip().split("\t")
                        species_list.append(item[0])
                        num_list.append(item[1])
                        percent_list.append(item[2])
                update_dict = {
                    "species_list": ":".join(species_list),
                    "species_num": ":".join(num_list),
                    "species_rate": ":".join(percent_list)
                }
            if update_dict:
                query_dict = {"split_id": split_id, "specimen_id": specimen_id}
                self.db["sg_split_clean_qc"].update(query_dict, {"$set": update_dict})

    def update_qc_stat_specimen_id(self, split_id, json_dir, product_type):
        """
        更新质控统计specimen_id
        """
        split_id = self.check_objectid(split_id)
        self.check_exists(json_dir)
        for f in os.listdir(json_dir):
            lib = f.split(".json")[0].split("--")[0]
            sample_name = f.split(".json")[0].split("--")[1]
            spe_result = self.db["sg_split_specimen"].find_one({"split_id": split_id, "library_number": lib, "specimen_name": sample_name})
            update_dict = {"specimen_id": spe_result["_id"]}
            self.db["sg_split_raw_qc"].update({"split_id": split_id, "library_number": lib, "specimen_name": sample_name}, {"$set": update_dict})

    def export_upload_list(self, split_id, type, upload_list):
        w = open(upload_list, "wb")
        split_id = self.check_objectid(split_id)
        if type == "library":
            results = self.db["sg_split_library"].find({"split_id": split_id})
            for result in results:
                work_path = result["work_path"].split(";")
                path = result["path"].split(";")
                for i in range(len(work_path)):
                    w.write(work_path[i] + "\t" + path[i] + "\n")
        elif type == "sample_split":
            results = self.db["sg_split_specimen"].find({"split_id": split_id})
            for result in results:
                work_path = result["work_path"].split(";")
                if result["product_type"] == "meta":
                    s3_path = result["clean_path"].split(";")
                else:
                    s3_path = result["raw_path"].split(";")
                for i in range(len(work_path)):
                    w.write(work_path[i] + "\t" + s3_path[i] + "\n")
        elif type == "qc":
            results = self.db["sg_split_specimen"].find({"split_id": split_id})
            for result in results:
                if result["product_type"] == "meta":
                    continue
                clean_work_path = result["clean_work_path"].split(";")
                clean_path = result["clean_path"].split(";")
                for i in range(len(clean_work_path)):
                    w.write(clean_work_path[i] + "\t" + clean_path[i] + "\n")
        w.close()

    def update_sg_split_specimen_merge(self, merge_id, output_dir, s3_upload_dir, operation_type):
        """
        更新合并样本的信息
        """
        self.check_exists(output_dir)
        merge_id = self.check_objectid(merge_id)
        md5sum_file = os.path.join(output_dir,"md5sum.txt")
        md5sum_dict = {}
        if os.path.exists(md5sum_file):
            with open(md5sum_file,"r") as mf:
                while 1:
                    line = mf.readline()
                    if not line:
                        break
                    fd = line.rstrip().split("  ")
                    md5sum_dict[fd[1]] = fd[0]
        raw_path, raw75_path, clean_path, raw_bytes, raw75_bytes, clean_bytes, raw_md5sum, clean_md5sum, raw75_md5sum = [], [], [], [], [], [], [], [], []
        for f in os.listdir(output_dir):
            path = os.path.join(output_dir, f)
            s3_path = os.path.join(s3_upload_dir, f)
            if f.endswith("R1.raw.fastq.gz"):
                raw_md5sum.insert(0,md5sum_dict[f])
                raw_path.insert(0, s3_path)
                raw_bytes.insert(0, str(os.path.getsize(path)))
            elif f.endswith("R2.raw.fastq.gz"):
                raw_md5sum.append(md5sum_dict[f])
                raw_path.append(s3_path)
                raw_bytes.append(str(os.path.getsize(path)))
            elif f.endswith("clean.1.fastq.gz"):
                clean_md5sum.insert(0, md5sum_dict[f])
                clean_path.insert(0, s3_path)
                clean_bytes.insert(0, str(os.path.getsize(path)))
            elif f.endswith("clean.2.fastq.gz"):
                clean_md5sum.append(md5sum_dict[f])
                clean_path.append(s3_path)
                clean_bytes.append(str(os.path.getsize(path)))
            elif (".clean.fastq"):  # 多样性
                clean_md5sum.append(md5sum_dict[f])
                clean_path.append(s3_path)
                clean_bytes.append(str(os.path.getsize(path)))
            elif f.endswith("clean.1.fastq"):  # 多样性
                clean_md5sum.insert(0, md5sum_dict[f])
                clean_path.insert(0, s3_path)
                clean_bytes.insert(0, str(os.path.getsize(path)))
            elif f.endswith("clean.2.fastq"):  # 多样性
                clean_md5sum.append(md5sum_dict[f])
                clean_path.append(s3_path)
                clean_bytes.append(str(os.path.getsize(path)))
            elif f.endswith(".fasta.gz"):
                clean_md5sum.append(md5sum_dict[f])
                clean_path.append(s3_path)
                clean_bytes.append(str(os.path.getsize(path)))
            elif f.endswith(".fq.gz"):
                raw75_md5sum.append(md5sum_dict[f])
                raw75_path.append(s3_path)
                raw75_bytes.append(str(os.path.getsize(path)))
        if operation_type == "merge":
            update_dict = {
                "raw75_md5sum":";".join(raw75_md5sum),
                "raw_md5sum":";".join(raw_md5sum),
                "clean_md5sum":";".join(clean_md5sum),
                "raw_path": ";".join(raw_path),
                "raw75_path": ";".join(raw75_path),
                "clean_path": ";".join(clean_path),
                "raw_bytes": ";".join(raw_bytes),
                "raw75_bytes": ";".join(raw75_bytes),
                "clean_bytes": ";".join(clean_bytes),
                "updated_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            }
            self.db["sg_split_specimen_merge"].update({"_id": merge_id}, {"$set": update_dict})
        elif operation_type == "rename":
            result = self.db["sg_split_specimen_rename"].find_one({"_id": merge_id})
            update_dict = {
                "re_raw_md5sum":result["raw_md5sum"],
                "re_raw_path": result["raw_path"],
                "re_clean_md5sum":";".join(clean_md5sum),
                "re_clean_path": ";".join(clean_path),
                "re_raw_bytes": result["raw_bytes"],
                "re_clean_bytes": ";".join(clean_bytes),
                "updated_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            }
            self.db["sg_split_specimen_rename"].update({"_id": merge_id}, {"$set": update_dict})

    def add_sg_split_library_qc(self, split_id, lib_qc_path):
        """
        导入meta 文库质控统计主表
        """
        split_id = self.check_objectid(split_id)
        self.check_exists(lib_qc_path)
        data_list = []
        lib_list = []
        with open(lib_qc_path, "rb") as f:
            lines = f.readlines()
            for line in lines[1:]:
                item = line.strip().split("\t")
                insert_data = {
                    "split_id": split_id,
                    "library_number": item[0],
                    "rank": item[1],
                    "q20": float(item[2]),
                    "q30": float(item[3]),
                    "raw_pair": int(item[4]),
                    "chimeric": int(item[5]),
                    "chimeric_rate": float(item[6]),
                    "valid_pair": int(item[7]),
                    "valid_rate": float(item[8]),
                    "pair_trim": int(item[9]),
                    "trim_rate": float(item[10]),
                    "pair_merge": int(item[11]),
                    "merge_rate": float(item[12]),
                    "seq_split": int(item[13]),
                    "split_rate": float(item[14]),
                    "high_quality_rate": float(item[15]),
                }
                data_list.append(insert_data)
                lib_list.append(item[0])
        results = self.db["sg_split_library"].find({"split_id": split_id})
        for result in results:
            if result["library_number"] not in lib_list:
                insert_data = {
                    "split_id": split_id,
                    "library_number": result["library_number"]
                }
                data_list.append(insert_data)
        if len(data_list) > 0:
            self.db["sg_split_library_qc"].insert_many(data_list)

    def add_sg_meta_verify_barcode(self):
        """
        导入meta barcode错乱验证主表
        """
        main_table_name = "barcode_disorder_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        split_id = "5d148cacd7f39d099218b292"
        library_number = "MJ190326_39"
        lib_type = "Illumina多样性文库"
        params_json = {
            "split_id": split_id,
            "library_number": library_number,
            "lib_type": lib_type,
        }
        params = json.dumps(params_json, sort_keys=True, separators=(',', ':'))
        data_list = []
        insert_data = {
            "name": main_table_name,
            "status": "start",
            "task_sn": main_table_name,
            "split_id": split_id,
            "library_number": library_number,
            "params": params,
            "type": "barcode_disorder",
            "desc": "多样性结果验证-barcode错乱",
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        }
        data_list.append(insert_data)
        self.db["sg_meta_verify_barcode"].insert_many(data_list)

    def add_sg_meta_verify_primer(self):
        """
        导入meta primer错配验证主表
        """
        main_table_name = "primer_mismatch_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        split_id = "5d148cacd7f39d099218b292"
        library_number = "MJ190326_39"
        params_json = {
            "split_id": split_id,
            "library_number": library_number,
        }
        params = json.dumps(params_json, sort_keys=True, separators=(',', ':'))
        data_list = []
        insert_data = {
            "name": main_table_name,
            "status": "start",
            "task_sn": main_table_name,
            "split_id": split_id,
            "library_number": library_number,
            "params": params,
            "type": "primer_mismatch",
            "desc": "多样性结果验证-primer错配",
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        }
        data_list.append(insert_data)
        self.db["sg_meta_verify_primer"].insert_many(data_list)

    def add_sg_meta_verify_barcode_detail(self, verify_id, barcode_path):
        """
        导入meta barcode错乱验证细节表
        """
        verify_id = self.check_objectid(verify_id)
        self.check_exists(barcode_path)
        data_list = []
        with open(barcode_path, "rb") as f:
            lines = f.readlines()
            for line in lines[1:]:
                item = line.strip().split("\t")
                insert_data = {
                    "verify_id": verify_id,
                    "barcode": item[0],
                    "reads": int(item[1]),
                    "f_barcode": item[2],
                    "r_barcode": item[3],
                }
                data_list.append(insert_data)
        self.db["sg_meta_verify_barcode_detail"].insert_many(data_list)

    def add_sg_meta_verify_primer_detail(self, verify_id, info_path):
        """
        导入meta primer错配验证细节表
        """
        verify_id = self.check_objectid(verify_id)
        self.check_exists(info_path)
        data_list = []
        with open(info_path, "rb") as f:
            lines = f.readlines()
            for line in lines[1:]:
                item = line.strip().split("\t")
                insert_data = {
                    "verify_id": verify_id,
                    "specimen_name": item[0].split("--")[-1],
                    "primer": item[1],
                    "reads": int(item[2])
                }
                data_list.append(insert_data)
        self.db["sg_meta_verify_primer_detail"].insert_many(data_list)

    def delete_sg_split_clean_qc(self, split_id):
        self.bind_object.logger.info("开始清除冗余数据")
        split_id = self.check_objectid(split_id)
        results1 = self.db["sg_split_specimen"].find({"split_id": ObjectId(split_id)})
        count1 = results1.count()
        results2 = self.db["sg_split_clean_qc"].find({"split_id": ObjectId(split_id)})
        count2 = results2.count()
        if count1 != count2:
            delete_num = 0
            sp_list = []
            for result in sorted(results2,key = lambda e:str(e.__getitem__("_id")),reverse=True):
                if str(result["specimen_id"]) in sp_list:
                    self.db["sg_split_clean_qc"].delete_one({"_id": result["_id"]})
                    delete_num = delete_num + 1
                    continue
                sp_list.append(str(result["specimen_id"]))
            self.bind_object.logger.info("清除sg_split_raw_qc表{}条数据".format(delete_num))

    def delete_sg_split_library_qc(self, split_id):
        # self.bind_object.logger.info("开始清除冗余数据")
        split_id = self.check_objectid(split_id)
        results1 = self.db["sg_split_library"].find({"split_id": ObjectId(split_id)})
        count1 = results1.count()
        results2 = self.db["sg_split_library_qc"].find({"split_id": ObjectId(split_id)})
        count2 = results2.count()
        if count1 != count2:
            delete_num = 0
            sp_list = []
            for result in sorted(results2,key = lambda e:str(e.__getitem__("_id")),reverse=True):
                if str(result["library_number"]) in sp_list:
                    self.db["sg_split_library_qc"].delete_one({"_id": result["_id"]})
                    delete_num = delete_num + 1
                    continue
                sp_list.append(str(result["library_number"]))
            # self.bind_object.logger.info("清除sg_split_library_qc表{}条数据".format(delete_num))

    def update_sg_split_specimen_status(self, coll_id, coll, status, desc):
        """
        更新合并样本的的运行状态
        """
        coll_id = self.check_objectid(coll_id)
        update_dict = {
            "status": status,
            "desc": desc
        }
        self.db[coll].update({"_id": coll_id}, {"$set": update_dict})

    def friendly_size(self, size):
        """
        资源转换
        """
        gb = 1000 * 1000 * 1000.0
        mb = 1000 * 1000.0
        kb = 1000.0
        if size > gb:
            new_size = round(float(size) / gb, 4)
            return str(new_size) + "G"
        else:
            new_size = round(float(size) / mb, 4)
            return str(new_size) + "M"

if __name__ == "__main__":
    a = DatasplitNew(None)
    # split_id = "6110f2d998ea795db322ea22"
    # lane_path = "/mnt/ilustre/users/sanger-dev/tsanger/workspace/20210810/LibrarySplit_CF3-20210802PE300-Mruijin_20210810_145153/Bcl2fastq/output/Reports/html/H735MDRXY/all/all/all/lane.html"
    # lane_barcode_path = "/mnt/ilustre/users/sanger-dev/tsanger/workspace/20210810/LibrarySplit_CF3-20210731bNovaSP_20210810_093405/Bcl2fastq/output/Reports/html/H735MDRXY/all/all/all/laneBarcode.html"
    # a.add_flowcell_summary1(split_id, lane_path, lane_barcode_path)
    split_id = "614d610c37959a083f62be76"
    # fastq_dir = "/mnt/ilustre/users/sanger-dev/tsanger/workspace/20210811/LibrarySplit_CF6-20210731bNovaSP_20210811_113550/output/2"
    # sample_sheet = "/mnt/ilustre/users/sanger-dev/tsanger/workspace/20210811/LibrarySplit_CF6-20210731bNovaSP_20210811_113550/2.sample_sheet.csv"
    # s3_upload_dir = "s3nb1://datasplit/2021/20210731bNovaSP/CF6-20210731bNovaSP_20210811_113550/2"
    # a.update_lib_path(split_id, fastq_dir, sample_sheet, s3_upload_dir)
    # a.delete_sg_split_library_qc(split_id)
    a.update_raw_sample_path(split_id)

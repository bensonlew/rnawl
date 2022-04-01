#!/usr/bin/env python
# -*- coding: utf-8 -*-
# __author__ == zhouxuan

import os
import re
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
import numpy as np
import pandas as pd


class BacGenomeAssemModule(Module):
    """
    对数据进行解压,需要根据不同样本文库，分别放置路径
    """
    def __init__(self, work_id):
        super(BacGenomeAssemModule, self).__init__(work_id)
        options = [
            {"name": "raw_dir", "type": "infile", "format": "bacgenome.raw_dir2", 'required': True},
            {"name": "qc", "type": "bool", "default": True},
            {"name": "qc_tool", "type": "string", "default": "fastp", 'choose': ['fastp', 'old_mode']},
            {"name": "dir", "type": "outfile", "format": "bacgenome.simple_dir"}
        ]
        self.add_option(options)
        self.zip_info = {}
        self.pacbio_info = {}
        self.tools = []
        self.pe_insert = {}  # 没有lib_name时的insert_size
        self.mp_insert = {}
        self.read_length = {}
        self.file_list = ["qc_a", "qc_b", "qc_stat", "bam_file", "third_file"]  # 解压后(fastp软件不需解压)，配置的文件，给后续的流程使用

    def check_options(self):
        """
        检查参数
        """
        return True

    def run(self):
        """
        三代数据转格式去存储
        二代数据进行解压缩，如果是fastp流程则不解压
        生成整理后的模块list文件，并提供完整的文件路径
        :return:
        """
        super(BacGenomeAssemModule, self).run()
        self.get_list_info()
        if not self.option("qc"):
            # 不做qc则必须保证碱基统计时所有的数据为解压缩后的
            self.run_ungiz(False)
        if self.option("qc_tool") == "old_mode":
            self.run_ungiz(False)
        elif self.option("qc_tool") == "fastp":
            # self.run_ungiz(True)
            self.run_ungiz(False)

    def add_zip_info(self, sample, lib_name, path_list):
        # 三代数据不运行此方法
        self.logger.info(self.zip_info)
        self.logger.info(sample)
        self.logger.info(lib_name)
        self.logger.info(path_list)
        path_num = len(path_list)
        if path_num != 2:
            self.set_error("样本%s,文库%s,双端数据个数不正确:%s" % (sample, lib_name, path_num))
        if self.zip_info.has_key(sample):
            if self.zip_info[sample].has_key(lib_name):
                if self.zip_info[sample][lib_name].has_key("1"):
                    self.zip_info[sample][lib_name]["1"] += " %s" % path_list[0]
                    self.zip_info[sample][lib_name]["2"] += " %s" % path_list[1]
                else:
                    self.zip_info[sample][lib_name]["1"] = path_list[0]
                    self.zip_info[sample][lib_name]["2"] = path_list[1]
            else:
                self.zip_info[sample][lib_name] = {
                    "1": path_list[0],
                    "2": path_list[1]
                }
        else:
            self.logger.info("add sample zip info ")
            self.zip_info[sample] = {
                lib_name: {
                    "1": path_list[0],
                    "2": path_list[1]
                }
            }

    def add_nanopore_info(self, sample, path_list):
        if self.zip_info.has_key(sample):
            if self.zip_info[sample].has_key("Nanopore"):
                self.zip_info[sample]["Nanopore"]["1"] += " %s" % (" ".join(path_list))
            else:
                self.zip_info[sample]["Nanopore"] = {"1":  " ".join(path_list)}
        else:
            self.zip_info[sample] = {
                "Nanopore": {
                    "1": " ".join(path_list)
                }
            }

    def add_pacbio_info(self, sample, path_list):
        if self.pacbio_info.has_key(sample):
            self.set_error("样品%s有重复的Pacbio数据信息" % sample)
        else:
            self.pacbio_info[sample] = " ".join(path_list)

    def parse(self, df):
        this_sample = df["Sample Name"]
        lib_name = df["Library Name"]
        lib_type = df["Library"]
        paths = df["File Name"].split(",")
        df["Insert Size(bp)"] = str(df["Insert Size(bp)"]).split(".")[0]
        df["Reads Length(bp)"] = str(df["Reads Length(bp)"]).split(".")[0]
        rel_paths = [os.path.join(self.option("raw_dir").prop['path'], i) for i in paths]
        if self.option("raw_dir").samples[this_sample][0] == "complete":
            if lib_type == "PE":
                if lib_name == "-":
                    self.add_zip_info(this_sample, "PE_" + str(df["Insert Size(bp)"]), rel_paths)
                    self.read_length[str(this_sample) + "_PE_" + str(df["Insert Size(bp)"])] = df["Reads Length(bp)"]
                else:
                    self.pe_insert[this_sample] = df["Insert Size(bp)"]
                    self.add_zip_info(this_sample, "PE_" + lib_name, rel_paths)
                    self.read_length[str(this_sample) + "_PE_" + lib_name] = df["Reads Length(bp)"]
            elif lib_type == "MP":
                self.mp_insert[this_sample] = df["Insert Size(bp)"]
                self.add_zip_info(this_sample, "MP_" , rel_paths)
                self.read_length[str(this_sample) + "_MP"] = df["Reads Length(bp)"]
            elif lib_type == "Nanopore":
                self.add_nanopore_info(this_sample, rel_paths)
            else:
                self.add_pacbio_info(this_sample, rel_paths)
        else:
            '''
            现在都是一样的
            '''
            if lib_type == "PE":
                if lib_name == "-":
                    self.add_zip_info(this_sample, "PE_" + str(df["Insert Size(bp)"]), rel_paths)
                    self.read_length[str(this_sample) + "_PE_" + str(df["Insert Size(bp)"])] = df["Reads Length(bp)"]
                else:
                    self.pe_insert[this_sample] = df["Insert Size(bp)"]
                    self.add_zip_info(this_sample, "PE_" + lib_name, rel_paths)
                    self.read_length[str(this_sample) + "_PE_" + lib_name] = df["Reads Length(bp)"]
            elif lib_type == "MP":
                self.mp_insert[this_sample] = df["Insert Size(bp)"]
                self.add_zip_info(this_sample, "MP", rel_paths)
                self.read_length[str(this_sample) + "_MP"] = df["Reads Length(bp)"]
            elif lib_type == "Nanopore":
                self.add_nanopore_info(this_sample, rel_paths)
            else:
                self.add_pacbio_info(this_sample, rel_paths)

    def get_list_info(self):
        list_txt = os.path.join(self.option("raw_dir").prop['path'], "list.txt")
        data = pd.read_table(list_txt)
        data.columns = ['Sample Name', 'File Name', 'Insert Size(bp)', 'Reads Length(bp)', 'Genome Size(Mb)', 'Library', 'Library Name']
        # modified by ghd @ 20190731,防止用新增的表头时报错
        # data["Insert Size(bp)"] = data["Insert Size(bp)"].apply(pd.to_numeric, errors='ignore').astype(np.int64, errors="ignore")
        # data["Reads Length(bp)"] = data["Reads Length(bp)"].apply(pd.to_numeric, errors="ignore").astype(np.int64, errors="ignore")
        # self.logger.info("sam1 debug")
        # self.logger.info(data)
        data.apply(self.parse, axis=1)

    def path_is_wrong(self, path_str):
        # 路径中的压缩是否是统一的？
        paths = path_str.split()
        if re.search('\.gz$', paths[0]) or re.search('\.gzip$', paths[0]):
            fisrt_gz_mark = "gz"
        else:
            fisrt_gz_mark = "ungz"
        for path in paths[1:]:
            if re.search('\.gz$', path) or re.search('\.gzip$', path):
                gz_mark = "gz"
            else:
                gz_mark = "ungz"
            if gz_mark != fisrt_gz_mark:
                return True

    def run_ungiz(self, zipornot):
        result_path = os.path.join(self.work_dir, "ungiz_dir")
        if not os.path.exists(result_path):
            os.mkdir(result_path)
        for sample in self.zip_info:
            for lib in self.zip_info[sample]:
                for d in self.zip_info[sample][lib]:
                    if self.path_is_wrong(self.zip_info[sample][lib][d]):
                        """
                        需要检查同一个tool里是否都是统一压缩形式
                        不是的话就跳过，整理得到的样品信息则不会出现此条数据
                        """
                        continue
                    gunzip_fastq = self.add_tool('bacgenome.fastq_ungz')
                    self.logger.info("sample: %s lib: %s" % (sample, lib))
                    self.logger.info(self.zip_info[sample])
                    gunzip_fastq.set_options({
                        "fastq": self.zip_info[sample][lib][d],
                        "sample_name": str(sample),
                        "direction": d,
                        "lib_type": lib,
                        "result_path": result_path,
                        "nozip": zipornot
                    })
                    self.tools.append(gunzip_fastq)
        self.run_bam()

    def run_bam(self):
        result_path = os.path.join(self.work_dir, "ungiz_dir")
        for sample in self.pacbio_info:
            pacbio_tool = self.add_tool('bacgenome.pacbio_convert')
            pacbio_tool.set_options({
                "files": self.pacbio_info[sample],
                "sample_name": str(sample),
                "fq_result_path": result_path
            })
            self.tools.append(pacbio_tool)
        if len(self.tools) > 1:
            self.on_rely(self.tools, self.set_output)
        elif len(self.tools) == 0:
            self.tools[0].on('end', self.set_output)
        for tool in self.tools:
            tool.run()

    def set_output(self):
        """
        要整理所有数据，包括三代的数据
        :return:
        """
        if self.option("qc") and self.option("qc_tool") == "fastp":
            write_file = os.path.join(self.work_dir, "qc_a")
        elif self.option("qc") and self.option("qc_tool") == "old_mode":
            write_file = os.path.join(self.work_dir, "qc_b")
        else:
            write_file = os.path.join(self.work_dir, "qc_stat")
        ungiz_dir = os.listdir(os.path.join(self.work_dir, "ungiz_dir"))
        ungiz_dir.sort()
        with open(write_file, "w") as write_record:
            for file in ungiz_dir:
                if "_PE_" in file:
                    sample,tmp_str = file.split("_PE_")
                    lib_type = "PE"
                    tmp_str2 = tmp_str.split(".")[0]
                    if tmp_str2.upper() == tmp_str2.lower():
                        insert_size = tmp_str2
                        lib_name = "-"
                    else:
                        if ".1.fq" in tmp_str:
                            lib_name = tmp_str.rstrip(".1.fq")
                        else:
                            lib_name = tmp_str.rstrip(".2.fq")
                        insert_size = self.pe_insert[sample]
                    read_len = self.read_length[sample + "_PE_" + tmp_str2]
                elif "_MP" in file:
                    sample,tmp_str = file.split("_MP")
                    lib_type = "MP"
                    insert_size = self.mp_insert[sample]
                    read_len = self.read_length[sample + "_MP"]
                    lib_name = "-"
                else:
                    continue
                if ".1.fq" in tmp_str:
                    write_record.write("%s\t%s\t%s\t%s\t%s\tl\t%s\n" % (sample, lib_type, lib_name, insert_size, read_len, file))
                else:
                    write_record.write("%s\t%s\t%s\t%s\t%s\tr\t%s\n" % (sample, lib_type, lib_name, insert_size, read_len, file))
        list_txt = os.path.join(self.option("raw_dir").prop['path'], "list.txt")
        size_dict = {}
        with open(list_txt, "r") as f:
            lines =f.readlines()
            for line in lines[1:]:
                lin = line.strip().split("\t")
                if lin[5] in ["Pacbio", "Nanopore"]:
                    size_dict[lin[0]] = lin[4]
        write_file = os.path.join(self.work_dir, "third_file")
        with open(write_file, "w") as write_record:
            for file in ungiz_dir:
                if "_Pacbio" in file:
                    sample,tmp_str = file.split("_Pacbio")
                    lib_type = "Pacbio"
                elif "_Nanopore" in file:
                    sample,tmp_str = file.split("_Nanopore")
                    lib_type = "Nanopore"
                else:
                    continue
                write_record.write("%s\t%s\t%s\t%s\n" % (sample, lib_type, file, size_dict[sample]))
        self.option("dir", os.path.join(self.work_dir, "ungiz_dir"))
        self.end()

    def end(self):
        super(BacGenomeAssemModule, self).end()

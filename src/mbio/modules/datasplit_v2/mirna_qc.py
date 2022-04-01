# -*- coding: utf-8 -*-
# __author__ = "zengjing"
# last_modify: 20171213

import os
import re
import time
import shutil
from collections import defaultdict
from biocluster.module import Module
from biocluster.core.exceptions import OptionError


class MirnaQcModule(Module):
    """
    miRNA质控的module,用于miRNA质控，多个fastq文件进行质控
    是否去除前三个碱基，去接头，去低值，去未知碱基序列、去过长过短序列
    文库+样本不能重复
    """
    def __init__(self, work_id):
        super(MirnaQcModule, self).__init__(work_id)
        options = [
            {"name": "list_file", "type": "infile", "format": "datasplit.list_file"},  # fastq文件及对应的样本
            {"name": "length", "type": "int", "default": 18},  # 最小序列长度，丢弃比此值更短的序列
            {"name": "adapter", "type": "string", "default": "TGGAATTCTCGGGTGCCAAGG"},  # 接头序列，如果微量建库，则改为AAAAAA
            {"name": "phred_score", "type": "string", "default": "20"},  # Phred得分（在0到40之间）
            {"name": "minlen", "type": "string", "default": "18"},  # 最短序列
            {"name": "max_length", "type": "string", "default": "32"},  # 最长序列
            {"name": "cut_left", "type":"bool", "default": False}  # miRNA是否要切除前3bp，保留前51bp
        ]
        self.add_option(options)
        self.end_times = 0

    def check_options(self):
        if not self.option("list_file").is_set:
            raise OptionError("必须设置list文件")

    def get_sample_info(self):
        self.sample_adapter = {}
        with open(self.option("list_file").prop["path"], "r") as f:
            for line in f:
                item = line.strip().split()
                if len(item) > 2:
                    self.sample_adapter[item[1]] = {"lib_type": item[2]}
                if len(item) == 4:
                    self.sample_adapter[item[1]]["adapter"] = item[4]

    def run_cut_left(self):
        for s in self.fastqs.keys():
            options = {
                "fastq": self.fastqs[s][0].split(";")[0],
                "cut_three": self.option("cut_left")
            }
            self.cut_left = self.add_tool("datasplit_v2.cut_left51")
            self.cut_left.set_options(options)
            self.cut_left.on("end", self.run_fastx_clipper, s)
            self.cut_left.on("end", self.run_fastq_stat, s)
            self.cut_left.on("end", self.set_output, "cut_left:" + s)
            self.cut_left.run()

    def run_fastq_stat(self, event):
        """
        质控前样本信息统计
        """
        options = {
            "fastq": event["bind_object"].option("cut_fastq").prop["path"]
        }
        self.fastq_stat = self.add_tool("datasplit_v2.fastq_stat")
        self.fastq_stat.set_options(options)
        self.fastq_stat.on("end", self.set_output, "fastq_stat:" + event["data"])
        self.fastq_stat.run()

    def run_fastx_clipper(self, event):
        options = {
            "length": self.option("length"),
            "adapter": self.option("adapter"),
            "fastq_s": event["bind_object"].option("cut_fastq").prop["path"]
        }
        if self.sample_adapter[event["data"]]:
            lib_type = self.sample_adapter[event["data"]]["lib_type"]
            if re.search("SR常量", lib_type):  # SR常量文库接头NNNNTGGAATTCTCGGGTGCCAAGG ，SR微量文库AAAAAAAAAAAAAAA
                options["adapter"] = "NNNNTGGAATTCTCGGGTGCCAAGG"
            elif re.search("SR微量", lib_type):
                options["adapter"] = "AAAAAAAAAAAAAAA"
            if "adapter" in self.sample_adapter[event["data"]].keys() and self.sample_adapter[event["data"]]["adapter"]:
                options["adapter"] = self.sample_adapter[event["data"]]["adapter"]
        self.fastx_clipper = self.add_tool("datasplit_v2.fastx_clipper")
        self.fastx_clipper.set_options(options)
        self.fastx_clipper.on("end", self.run_dynamic_trim, event["data"])
        self.fastx_clipper.on("end", self.set_output, "fastx_clipper:" + event["data"])
        self.fastx_clipper.run()

    def run_dynamic_trim(self, event):
        options = {
            "fastq": event["bind_object"].option("clip_s"),
            "phred_score": self.option("phred_score")
        }
        self.dynamic_trim = self.add_tool("datasplit_v2.dynamic_trim")
        self.dynamic_trim.set_options(options)
        self.dynamic_trim.on("end", self.run_fastq_fasta, event["data"])
        self.dynamic_trim.run()

    def run_fastq_fasta(self, event):
        self.logger.info(event["bind_object"].option("out_fastq").prop["path"])
        options = {
           "fastq": event["bind_object"].option("out_fastq").prop["path"]
        }
        self.fastq_fasta = self.add_tool("datasplit_v2.fastq_to_fasta")
        self.fastq_fasta.set_options(options)
        self.fastq_fasta.on("end", self.set_output, "fastq_to_fasta:" + event["data"])
        self.fastq_fasta.on("end", self.run_fasta_length, event["data"])
        self.fastq_fasta.run()

    def run_fasta_length(self, event):
        """
        质控后统计，计算fasta的长度分布情况
        """
        options = {
            "fasta": os.path.join(event["bind_object"].output_dir, "fasta")
        }
        self.fasta_length = self.add_tool("datasplit_v2.fasta_length")
        self.fasta_length.set_options(options)
        self.fasta_length.on("end", self.set_output, "fasta_length:" + event["data"])
        self.fasta_length.run()

    def link(self, old, new):
        if os.path.exists(new):
            os.remove(new)
        os.link(old, new)

    def set_output(self, event):
        """
        输出文件：
        每个样本fastq_stat的结果写到fastq_stat.xls
        每个样本18nt、32nt质控的统计结果写到Sample_QC_stat.xls
        """
        obj = event["bind_object"]
        if event["data"].startswith("cut_left"):
            s = event["data"].split("cut_left:")[1]
            for f in os.listdir(obj.output_dir):
                if f.endswith(".gz"):
                    fastq_dir = os.path.join(self.output_dir, "fastq")
                    if not os.path.exists(fastq_dir):
                        os.mkdir(fastq_dir)
                    old = os.path.join(obj.output_dir, f)
                    new = os.path.join(fastq_dir, s + ".fq.gz")
                    self.link(old, new)
        elif event["data"].startswith("fastq_stat"):
            s = event["data"].split("fastq_stat:")[1]
            with open(os.path.join(obj.output_dir, "fastq_stat.xls"), "r") as f:
                lines = f.readlines()
                self.raw_stat[s] = lines[1].strip().split("\t")
        elif event["data"].startswith("fastx_clipper"):
            s = event["data"].split("fastx_clipper:")[1]
            old = os.path.join(obj.output_dir, "Sample_QC_stat.xls")
            with open(old, "r") as f:
                lines = f.readlines()
                self.qc_stat[s]["cliper"] =lines[1].strip().split("\t")
        elif event["data"].startswith("fastq_to_fasta"):
            s = event["data"].split("fastq_to_fasta:")[1]
            old = os.path.join(obj.output_dir, "reads_stat.xls")
            with open(old, "r") as f:
                lines = f.readlines()
                self.qc_stat[s]["fasta"] = lines[1].strip().split("\t")
            old = os.path.join(obj.output_dir, "fasta.gz")
            fasta_dir = os.path.join(self.output_dir, "fasta")
            if not os.path.exists(fasta_dir):
                os.mkdir(fasta_dir)
            new = os.path.join(fasta_dir, s + ".fasta.gz")
            self.link(old, new)
        elif event["data"].startswith("fasta_length"):
            fasta_length_dir = os.path.join(self.output_dir, "fasta_length")
            if not os.path.exists(fasta_length_dir):
                os.mkdir(fasta_length_dir)
            fasta_length = os.path.join(fasta_length_dir, event["data"].split("fasta_length:")[1] + ".fasta_length.xls")
            self.link(os.path.join(obj.output_dir, "fasta_length.xls"), fasta_length)
        self.end_times += 1
        if self.end_times == 5 * len(self.fastqs):
            with open(self.output_dir + "/Sample_QC_stat.xls", "w") as w:
                w.write("Sample\tRaw_reads\tAdapter_only\tN_reads\t<18nt\t>32nt\tClean_reads\tAdapter%\n")
                for s in self.qc_stat.keys():
                    # w.write(s + "\t" + "\t".join(self.qc_stat[s]["cliper"][1:4]) + "\t" + "\t".join(self.qc_stat[s]["fasta"]))
                    # w.write("\t" + self.qc_stat[s]["cliper"][-1] + "\n")
                    w.write(s + "\t" + "\t".join(self.qc_stat[s]["cliper"][1:4]) + "\t")
                    w.write(str(int(self.qc_stat[s]["fasta"][0]) + int(self.qc_stat[s]["cliper"][4])) + "\t")
                    w.write(str(int(self.qc_stat[s]["fasta"][1]) + int(self.qc_stat[s]["cliper"][5])) + "\t")
                    w.write(self.qc_stat[s]["fasta"][2] + "\t" + self.qc_stat[s]["cliper"][-1] + "\n")
            with open(self.output_dir + "/fastq_stat.xls", "w") as w:
                w.write("#Sample_ID\tTotal_Reads\tTotal_Bases\tTotal_Reads_with_Ns\tN_Reads%\tA%\tT%\tC%\tG%\tN%\tError%\tQ20%\tQ30%\tGC%\n")
                for s in self.raw_stat.keys():
                    w.write(s + "\t" + "\t".join(self.raw_stat[s][1:]) + "\n")
            self.end()

    def run(self):
        super(MirnaQcModule, self).run()
        self.fastqs = self.option("list_file").prop["samples"]
        self.get_sample_info()
        self.qc_stat = defaultdict(dict)
        self.raw_stat = defaultdict(list)
        self.run_cut_left()

    def end(self):
        super(MirnaQcModule, self).end()

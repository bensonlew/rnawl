# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# last_modified: 20171130

import os
import re
from biocluster.config import Config
from biocluster.module import Module
from biocluster.core.exceptions import OptionError


class SingleMirnaQcModule(Module):
    """
    miRNA质控的module,用于miRNA质控,单个fastq文件进行质控
    是否去除前三个碱基，去接头，去低值，去未知碱基序列、去过长过短序列
    """
    def __init__(self, work_id):
        super(SingleMirnaQcModule, self).__init__(work_id)
        options = [
            {"name": "fastq", "type": "infile", "format": "sequence.fastq"},  # fastq文件
            {"name": "length", "type": "int", "default": 18},  # 最小序列长度，丢弃比此值更短的序列
            {"name": "adapter", "type": "string", "default": 'TGGAATTCTCGGGTGCCAAGG'},  # 接头序列，如果微量建库，则改为AAAAAA
            {"name": "phred_score", "type": "string", "default": "20"},  # Phred得分（在0到40之间）
            {"name": "min_length", "type": "string", "default": "18"},  # 最短序列
            {"name": "max_length", "type": "string", "default": "32"},  # 最长序列
            {"name": "cut_left", "type": "string", "default": "True"}  # miRNA是否要切除前3bp，保留前51bp
        ]
        self.add_option(options)
        self.end_times = 0

    def check_options(self):
        if not self.option("fastq").is_set:
            raise OptionError("必须设置fastq文件")

    def run_cut_left(self):
        options = {
            "fastq": self.option("fastq"),
            "cut_three": self.option("cut_left")
        }
        self.cut_left = self.add_tool("datasplit.cut_left51")
        self.cut_left.set_options(options)
        self.cut_left.on("end", self.run_fastx_clipper)
        self.cut_left.on("end", self.run_fastq_stat)
        self.cut_left.on("end", self.run_fastx)
        self.cut_left.run()

    def run_fastq_stat(self):
        options = {
            "fastq": self.cut_left.option("cut_fastq").prop["path"]
        }
        self.fastq_stat = self.add_tool("datasplit.fastq_stat")
        self.fastq_stat.set_options(options)
        self.fastq_stat.on('end', self.set_output, 'fastq_stat')
        self.fastq_stat.run()

    def run_fastx(self):
        options = {
            "fastq": self.cut_left.option("cut_fastq").prop["path"]
        }
        self.fastx = self.add_tool("datasplit.fastx_v2")
        self.fastx.set_options(options)
        self.fastx.on('end', self.set_output, 'fastx')
        self.fastx.run()

    def run_fastx_clipper(self):
        options = {
            "length": self.option("length"),
            "adapter": self.option("adapter")
        }
        options["fastq_s"] = self.cut_left.option("cut_fastq")
        self.fastx_clipper = self.add_tool("datasplit.fastx_clipper")
        self.fastx_clipper.set_options(options)
        self.fastx_clipper.on("end", self.run_dynamic_trim)
        self.fastx_clipper.run()

    def run_dynamic_trim(self):
        options = {
            "fastq": self.fastx_clipper.option("clip_s"),
            "phred_score": self.option("phred_score")
        }
        self.dynamic_trim = self.add_tool("datasplit.dynamic_trim")
        self.dynamic_trim.set_options(options)
        self.dynamic_trim.on("end", self.run_fastq_fasta)
        self.dynamic_trim.run()

    def run_fastq_fasta(self):
        options = {
           "fastq": self.dynamic_trim.option("out_fastq")
        }
        self.fastq_fasta = self.add_tool("datasplit.fastq_to_fasta")
        self.fastq_fasta.set_options(options)
        self.fastq_fasta.on('end', self.set_output, 'fastq_to_fasta')
        # self.fastq_fasta.on("end", self.run_trim_seq)
        self.fastq_fasta.run()

    # def run_trim_seq(self):
    #     options = {
    #         "fasta": os.path.join(self.fastq_fasta.output_dir, "fasta"),
    #         "min_length": self.option("min_length"),
    #         "max_length": self.option("max_length")
    #     }
    #     self.trim_seq = self.add_tool("datasplit.trim_seq")
    #     self.trim_seq.set_options(options)
    #     self.trim_seq.on('end', self.set_output, 'trim_seq')
    #     self.trim_seq.run()

    def set_output(self, event):
        obj = event["bind_object"]
        if event["data"] == "fastq_to_fasta":
            stat1 = self.fastx_clipper.output_dir + "/Sample_QC_stat.xls"
            stat2 = self.fastq_fasta.output_dir + "/reads_stat.xls"
            with open(stat1, "r") as f1, open(stat2, "r") as f2, open(self.output_dir + "/Sample_QC_stat.xls", "w") as w:
                lines1 = f1.readlines()
                lines2 = f2.readlines()
                item1 = lines1[1].strip().split("\t")
                item2 = lines2[1].strip().split("\t")
                less_reads = int(item1[4]) + int(item2[0])
                w.write("Sample\tRaw_reads\tAdapter_only\tN_reads\t<18nt\t>32nt\tClean_reads\tAdapter%\n")
                w.write(item1[0] + "\t" + item1[1] + "\t" + item1[2] + "\t" + item1[3] + "\t" + str(less_reads) + "\t")
                w.write(item2[1] + "\t" + item2[2] + "\t" + item1[-1] + "\n")
            if os.path.exists(os.path.join(self.output_dir, "clean.fasta.gz")):
                os.remove(os.path.join(self.output_dir, "clean.fasta.gz"))
            os.link(os.path.join(self.fastq_fasta.output_dir, "fasta.gz"), os.path.join(self.output_dir, "clean.fasta.gz"))
        else:
            for f in os.listdir(obj.output_dir):
                if os.path.exists(os.path.join(self.output_dir, f)):
                    os.remove(os.path.join(self.output_dir, f))
                os.link(os.path.join(obj.output_dir, f), os.path.join(self.output_dir, f))
        self.end_times += 1
        if self.end_times == 3:
            self.end()

    def run(self):
        super(SingleMirnaQcModule, self).run()
        self.run_cut_left()

    def end(self):
        super(SingleMirnaQcModule, self).end()

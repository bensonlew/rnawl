# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# last_modify: 20171213

import os
import re
import time
import shutil
from biocluster.module import Module
from biocluster.core.exceptions import OptionError


class MirnaQcModule(Module):
    """
    miRNA质控的module,用于miRNA质控，多个fastq文件进行质控
    是否去除前三个碱基，去接头，去低值，去未知碱基序列、去过长过短序列
    """
    def __init__(self, work_id):
        super(MirnaQcModule, self).__init__(work_id)
        options = [
            {"name": "list_file", "type": "infile", "format": "datasplit.list_file"},  # fastq文件及对应的样本
            {"name": "length", "type": "int", "default": 18},  # 最小序列长度，丢弃比此值更短的序列
            {"name": "adapter", "type": "string", "default": 'TGGAATTCTCGGGTGCCAAGG'},  # 接头序列，如果微量建库，则改为AAAAAA
            {"name": "phred_score", "type": "string", "default": "20"},  # Phred得分（在0到40之间）
            {"name": "minlen", "type": "string", "default": "18"},  # 最短序列
            {"name": "max_length", "type": "string", "default": "32"},  # 最长序列
            {"name": "cut_left", "type": "string", "default": "True"}  # miRNA是否要切除前3bp，保留前51bp
        ]
        self.add_option(options)
        self.step.add_steps("single_mirna_qc")
        self.fastqs, self.cut_fastq = {}, {}
        self.end_times = 0

    def check_options(self):
        if not self.option("list_file").is_set:
            raise OptionError("必须设置list文件")

    def run_mirna_qc(self):
        options = {
            "length": self.option("length"),
            "adapter": self.option("adapter"),
            "phred_score": self.option("phred_score"),
            "min_length": self.option("minlen"),
            "max_length": self.option("max_length"),
            "cut_left": self.option("cut_left")
        }
        modules = []
        for s in self.fastqs.keys():
            options["fastq"] = self.fastqs[s][0].split(";")[0]
            self.single_mirna_qc = self.add_module("datasplit.single_mirna_qc")
            self.single_mirna_qc.set_options(options)
            self.single_mirna_qc.on('end', self.set_output, 'mirna_qc_{}'.format(s))
            modules.append(self.single_mirna_qc)
        if len(modules) == 1:
            modules[0].on("end", self.end)
        else:
            self.on_rely(modules, self.end)
        for m in modules:
            m.run()

    def link_file(self, old, new):
        if os.path.exists(new):
            os.remove(new)
        os.link(old, new)

    def set_output(self, event):
        obj = event["bind_object"]
        olddir = obj.output_dir
        sample = event["data"].split("mirna_qc_")[1]
        for f in os.listdir(olddir):
            if re.search(r".+.fasta.gz", f):
                f1 = sample + "_clean.fasta.gz"
                self.link_file(os.path.join(olddir, f), os.path.join(self.output_dir, "fasta/" + f1))
            elif re.search(r"fastq_stat.xls", f):
                with open(os.path.join(olddir, f), "r") as f, open(self.fastq_stat, "a+") as w:
                    lines = f.readlines()
                    for line in lines[1:]:
                        item = line.strip().split("\t")[1:]
                        w.write(sample + "\t" + "\t".join(item) + "\n")
            elif re.search(r".+.xls", f):
                f1 = sample + "_qc_stat.xls"
                self.link_file(os.path.join(olddir, f), os.path.join(self.output_dir, "fastq_stat/fastq_stat/" + f1))
            else:
                if f.endswith("fastxstat"):
                    f1 = sample + "_fastxstat"
                elif f.endswith("q20q30"):
                    f1 = sample + "_q20q30"
                else:
                    f1 = sample + ":" + f
                self.link_file(os.path.join(olddir, f), os.path.join(self.output_dir, "fastq_stat/fastx/" + f1))

    def run(self):
        super(MirnaQcModule, self).run()
        self.fastqs = self.option("list_file").prop["samples"]
        for f in os.listdir(self.output_dir):
            shutil.rmtree(os.path.join(self.output_dir, f))
        os.mkdir(os.path.join(self.output_dir, "fastq_stat"))
        os.mkdir(os.path.join(self.output_dir, "fastq_stat/fastq_stat"))
        os.mkdir(os.path.join(self.output_dir, "fastq_stat/fastx"))
        os.mkdir(os.path.join(self.output_dir, "fasta"))
        self.fastq_stat = os.path.join(self.output_dir, "fastq_stat/fastq_stat/fastq_stat.xls")
        with open(self.fastq_stat, "w") as w:
            w.write("#Sample_ID\tTotal_Reads\tTotal_Bases\tTotal_Reads_with_Ns\tN_Reads%\tA%\tT%\tC%\tG%\tN%\tError%\tQ20%\tQ30%\tGC%\n")
        self.run_mirna_qc()

    def end(self):
        super(MirnaQcModule, self).end()

# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# last_modify: 20171214

import os
import re
import shutil
from biocluster.module import Module
from biocluster.core.exceptions import OptionError


class CleanDataStatModule(Module):
    """
    质控后数据统计，统计冗余度，fastq_dup.py(输入为文件夹时文件夹内需要有list.txt)
    统计碱基质量信息及Q20、Q30，fastx_quality_stats、q20q30_stat.py(输入为文件夹时文件夹内需要有list.txt)
    统计fastq序列基本信息，FastqStat.jar(输入为文件夹时文件夹内需要有list.txt)
    DNA/常规RNA时加：随机抽取10000条序列，比对到NT数据库，进行统计
    常规RNA时加：随机抽取10000条序列，比对到Rfam数据库，进行统计
    miRNA时只做：统计总的fasta reads数，序列长度统计(质控时fastx_clipper的输出统计文件)
    """
    def __init__(self, work_id):
        super(CleanDataStatModule, self).__init__(work_id)
        options = [
            {"name": "list_file", "type": "infile", "format": "datasplit.list_file"},  # 要进行质控的文件的list,第一列文件路径，二列样本，三列r/l
            {"name": "project_type", "type": "string"},  # list_file里的文件对应的项目类型
        ]
        self.add_option(options)
        self.step.add_steps("dup", "stat", "fastx", "fastq_blast_nt", "fastq_blast_rfam", "fasta_length")
        self.samples = {}
        self.start_times = 0
        self.end_times = 0

    def set_step(self, event):
        if "start" in event["data"].keys():
            event["data"]["start"].start()
        if "end" in event["data"].keys():
            event["data"]["end"].finish()
        self.step.update()

    def check_options(self):
        if not self.option("list_file").is_set:
            raise OptionError("必须设置进行质控的文件list")
        self.samples = self.option("list_file").prop["samples"]
        if not self.option("project_type"):
            raise OptionError("必须设置所属的项目类型")
        if self.option("project_type") not in ["rna", "mirna", "ncrna", "dna", "meta", "meta_genomic", "microbial_genome"]:
            raise OptionError("项目类型{}必须为rna/mirna/ncrna/dna/meta/meta_genomic/microbial_genome".format(self.option("project_type")))

    def run_fastq_dup(self):
        for s in self.samples.keys():
            if len(self.samples[s]) == 1:
                options = {
                    "fastq_s": self.samples[s][0],
                    "fq_type": "SE"
                }
            else:
                options = {
                    "fastq_l": self.samples[s][0],
                    "fastq_r": self.samples[s][1],
                    "fq_type": "PE"
                }
            self.fastq_dup = self.add_tool("datasplit.fastq_dup")
            self.fastq_dup.set_options(options)
            self.fastq_dup.on("start", self.set_step, {"start": self.step.dup})
            self.fastq_dup.on("end", self.set_step, {"end": self.step.dup})
            self.fastq_dup.on('end', self.set_output, 'fastq_dup_{}'.format(s))
            self.fastq_dup.run()
            self.start_times += 1

    def run_fastq_stat(self):
        options = {
            "list_file": self.option("list_file")
        }
        self.fastq_stat = self.add_tool("datasplit.fastq_stat")
        self.fastq_stat.set_options(options)
        self.fastq_stat.on("start", self.set_step, {"start": self.step.stat})
        self.fastq_stat.on("end", self.set_step, {"end": self.step.stat})
        self.fastq_stat.on('end', self.set_output, 'fastq_stat')
        self.fastq_stat.run()
        self.start_times += 1

    def run_fastx(self):
        for s in self.samples.keys():
            if len(self.samples[s]) == 2:
                options = {
                    "fastq": self.samples[s][1]
                }
                self.fastx = self.add_tool("datasplit.fastx_v2")
                self.fastx.set_options(options)
                self.fastx.on("start", self.set_step, {"start": self.step.fastx})
                self.fastx.on("end", self.set_step, {"end": self.step.fastx})
                self.fastx.on('end', self.set_output, 'fastx_{}_r'.format(s))
                self.fastx.run()
                self.start_times += 1
                options = {
                    "fastq": self.samples[s][0]
                }
                self.fastx = self.add_tool("datasplit.fastx_v2")
                self.fastx.set_options(options)
                self.fastx.on("start", self.set_step, {"start": self.step.fastx})
                self.fastx.on("end", self.set_step, {"end": self.step.fastx})
                self.fastx.on('end', self.set_output, 'fastx_{}_l'.format(s))
                self.fastx.run()
                self.start_times += 1
            else:
                options = {
                    "fastq": self.samples[s][0]
                }
                self.fastx = self.add_tool("datasplit.fastx_v2")
                self.fastx.set_options(options)
                self.fastx.on("start", self.set_step, {"start": self.step.fastx})
                self.fastx.on("end", self.set_step, {"end": self.step.fastx})
                self.fastx.on('end', self.set_output, 'fastx_{}'.format(s))
                self.fastx.run()
                self.start_times += 1

    def run_nt_stat(self):
        for s in self.samples.keys():
            if len(self.samples[s]) == 2:
                options = {
                    "fastq_l": self.samples[s][0],
                    "fastq_r": self.samples[s][1]
                }
            else:
                options = {
                    "fastq": self.samples[s][0]
                }
            options["database"] = "nt"
            self.fastq_blast_nt = self.add_module("datasplit.single_fastq_blast")
            self.fastq_blast_nt.set_options(options)
            self.fastq_blast_nt.on("start", self.set_step, {"start": self.step.fastq_blast_nt})
            self.fastq_blast_nt.on("end", self.set_step, {"end": self.step.fastq_blast_nt})
            self.fastq_blast_nt.on('end', self.set_output, 'fastq_blast_nt_{}'.format(s))
            self.fastq_blast_nt.run()
            self.start_times += 1

    def run_rfam_stat(self):
        for s in self.samples.keys():
            if len(self.samples[s]) == 2:
                options = {
                    "fastq_l": self.samples[s][0],
                    "fastq_r": self.samples[s][1]
                }
            else:
                options = {
                    "fastq": self.samples[s][0]
                }
            options["database"] = "rfam"
            self.fastq_blast_rfam = self.add_module("datasplit.single_fastq_blast")
            self.fastq_blast_rfam.set_options(options)
            self.fastq_blast_rfam.on("start", self.set_step, {"start": self.step.fastq_blast_rfam})
            self.fastq_blast_rfam.on("end", self.set_step, {"end": self.step.fastq_blast_rfam})
            self.fastq_blast_rfam.on('end', self.set_output, 'fastq_blast_rfam_{}'.format(s))
            self.fastq_blast_rfam.run()
            self.start_times += 1

    def run_fasta_length(self):
        for s in self.samples.keys():
            if len(self.samples[s]) == 2:
                self.set_error("miRNA有双端？？？")
            else:
                options = {
                    "fasta": self.samples[s][0]
                }
            self.fasta_length = self.add_tool("datasplit.fasta_length")
            self.fasta_length.set_options(options)
            self.fasta_length.on("start", self.set_step, {"start": self.step.fasta_length})
            self.fasta_length.on("end", self.set_step, {"end": self.step.fasta_length})
            self.fasta_length.on('end', self.set_output, 'fasta_length_{}'.format(s))
            self.fasta_length.run()
            self.start_times += 1

    def set_output(self, event):
        self.logger.info("开始set output")
        obj = event["bind_object"]
        self.end_times += 1
        m1 = re.match(r"fastq_dup_(.+)", event["data"])
        m2 = re.match(r"fastx_(.+)", event["data"])
        m3 = re.match(r"fastq_blast_nt_(.+)", event["data"])
        m4 = re.match(r"fastq_blast_rfam_(.+)", event["data"])
        m5 = re.match(r"fasta_length_(.+)", event["data"])
        if m1:
            new_dir = os.path.join(self.output_dir, "fastq_dup")
            if not os.path.exists(new_dir):
                os.mkdir(new_dir)
            for f in os.listdir(obj.output_dir):
                if os.path.exists(new_dir + "/" + m1.group(1) + "_dup.xls"):
                    os.remove(new_dir + "/" + m1.group(1) + "_dup.xls")
                os.link(obj.output_dir + "/" + f, new_dir + "/" + m1.group(1) + "_dup.xls")
        if m2:
            new_dir = os.path.join(self.output_dir, "fastx")
            if not os.path.exists(new_dir):
                os.mkdir(new_dir)
            for f in os.listdir(obj.output_dir):
                n1 = re.match(r".+fastxstat$", f)
                n2 = re.match(r".+q20q30$", f)
                if n1:
                    if os.path.exists(new_dir + "/" + m2.group(1) + "_fastxstat"):
                        os.remove(new_dir + "/" + m2.group(1) + "_fastxstat")
                    os.link(obj.output_dir + "/" + f, new_dir + "/" + m2.group(1) + "_fastxstat")
                if n2:
                    if os.path.exists(new_dir + "/" + m2.group(1) + "_q20q30"):
                        os.remove(new_dir + "/" + m2.group(1) + "_q20q30")
                    os.link(obj.output_dir + "/" + f, new_dir + "/" + m2.group(1) + "_q20q30")
        if event["data"] == "fastq_stat":
            new_dir = os.path.join(self.output_dir, "fastq_stat")
            if os.path.exists(new_dir):
                shutil.rmtree(new_dir)
            os.mkdir(new_dir)
            for f in os.listdir(obj.output_dir):
                if os.path.exists(new_dir + "/" + f):
                    os.remove(new_dir + "/" + f)
                os.link(obj.output_dir + "/" + f, new_dir + "/" + f)
        if m3:
            new_dir = os.path.join(self.output_dir, "blast_nt")
            if not os.path.exists(new_dir):
                os.mkdir(new_dir)
            if os.path.exists(new_dir + "/" + m3.group(1) + "_nt_species.xls"):
                os.remove(new_dir + "/" + m3.group(1) + "_nt_species.xls")
            os.link(obj.output_dir + "/nt_species_stat.xls", new_dir + "/" + m3.group(1) + "_nt_species.xls")
        if m4:
            new_dir = os.path.join(self.output_dir, "blast_rfam")
            if not os.path.exists(new_dir):
                os.mkdir(new_dir)
            if os.path.exists(new_dir + "/" + m4.group(1) + "_rfam_summary.xls"):
                os.remove(new_dir + "/" + m4.group(1) + "_rfam_summary.xls")
            os.link(obj.output_dir + "/rfam_summary.xls", new_dir + "/" + m4.group(1) + "_rfam_summary.xls")
        if m5:
            new_dir = os.path.join(self.output_dir, "fasta_length")
            if not os.path.exists(new_dir):
                os.mkdir(new_dir)
            for f in os.listdir(obj.output_dir):
                if os.path.exists(new_dir + "/" + m5.group(1) + "_length.xls"):
                    os.remove(new_dir + "/" + m5.group(1) + "_length.xls")
                os.link(obj.output_dir + "/" + f, new_dir + "/" + m5.group(1) + "_length.xls")
        if self.start_times == self.end_times:
            self.logger.info("结束set output")
            self.end()

    def run(self):
        super(CleanDataStatModule, self).run()
        if self.option("project_type") != "mirna":
            self.run_fastq_dup()
            self.run_fastq_stat()
            self.run_fastx()
        if self.option("project_type") == "rna":
            self.run_nt_stat()
            self.run_rfam_stat()
        if self.option("project_type") == "dna":
            self.run_nt_stat()
        if self.option("project_type") == "mirna":
            self.run_fasta_length()

    def end(self):
        super(CleanDataStatModule, self).end()

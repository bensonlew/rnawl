# -*- coding: utf-8 -*-
# __author__ = 'shicaiping'
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import re
import ConfigParser
import unittest
from Bio import SeqIO
import pandas as pd

class SrnaStatAgent(Agent):
    """
    Rfam比对结果统计
    """
    def __init__(self, parent):
        super(SrnaStatAgent, self).__init__(parent)
        options = [
            {"name": "rfam_summary", "type": "string", "default": ""},  # rfam注释结果文件
            {"name": "repeat_summary", "type": "string", "default": ""},  # repeatmasker注释结果文件
            {"name": "exon_summary", "type": "string", "default": ""},  # exon注释结果文件
            {"name": "intron_summary", "type": "string", "default": ""},  # intron注释结果文件
            {"name": "known_mirna_count", "type": "string", "default": ""},  # 已知miRNA定量结果文件，为了获取全部样本名
            {"name": "novel_mirna_count", "type": "string", "default": ""},  # 新预测miRNA定量结果文件，为了获取全部样本名
            {"name": "known_mrd", "type": "string", "default": ""},  # 已知miRNA鉴定表
            {"name": "novel_mrd", "type": "string", "default": ""},  # 预测miRNA鉴定表
            {"name": "input_fa", "type": "infile", "format": "small_rna.fasta"}, # 质控后fasta文件
            {"name": "config", "type": "infile", "format": "small_rna.common"},  # 样本名转换文件
            {"name": "category", "type": "string", "default": ""},  # 物种分类，Animal or Plant
            {"name": "method", "type": "string", "default": ""},  # mirdeep2, mireap, mirdp2, mir_prefer
            {"name": "list", "type": "string", "default": ""},  # 样本列表，用来决定样本展示顺序
            {"name": "qc_output", "type": "string", "default": ""}, # 质控结果目录
            {"name": "input_type", "type": "string", "default": "raw"},  # srna模块分析的输入序列，raw为质控后，clean为mapping后
        ]
        self.add_option(options)
        self.step.add_steps("srna_stat")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.srna_stat.start()
        self.step.update()

    def stepfinish(self):
        self.step.srna_stat.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检测函数
        """
        if not os.path.exists(self.option("rfam_summary")):
            raise OptionError("rfam_summary文件不存在")
        if not os.path.exists(self.option("repeat_summary")):
            raise OptionError("repeat_summary文件不存在")
        if not os.path.exists(self.option("exon_summary")):
            raise OptionError("exon_summary文件不存在")
        if not os.path.exists(self.option("intron_summary")):
            raise OptionError("intron_summary文件不存在")
        if not os.path.exists(self.option("input_fa").prop["path"]):
            raise OptionError("input_fa文件不存在")
        if not os.path.exists(self.option("config").prop["path"]):
            raise OptionError("config文件不存在")
        return True

    def set_resource(self):
        """
        设置所需资源，需在类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 2
        self._memory = '10G'

    def end(self):
        super(SrnaStatAgent, self).end()

class SrnaStatTool(Tool):
    def __init__(self, config):
        super(SrnaStatTool, self).__init__(config)
        self.python = "miniconda2/bin/python"

    def run(self):
        super(SrnaStatTool, self).run()
        self.run_srna_stat()
        self.end()

    def parse_config(self, file=None, section=None, name=None):
        config = ConfigParser.ConfigParser()
        config.read(file)
        return config.get(section, name)

    def parse_config2(self, file=None):
        config = ConfigParser.ConfigParser()
        config.read(file)
        return config

    def run_srna_stat(self):
        self.logger.info("开始进行srna注释统计")
        self.samples = []
        with open(self.option("list"), "r") as f:
            for line in f:
                if line.strip().split("\t")[1].strip() not in self.samples:
                    self.samples.append(line.strip().split("\t")[1].strip())

        ## mirna stat
        dict_known = {}
        dict_novel = {}
        if self.option("known_mirna_count"):
            with open(self.option("known_mirna_count"), "r") as f1:
                known_mirna = {}
                head = f1.readline()
                samples = head.strip().split("\t")[1:]
                for sample in samples:
                    dict_known[sample] = 0
                for line in f1:
                    mirna_id = line.strip().split("\t")[0]
                    if mirna_id not in known_mirna:
                        known_mirna[mirna_id] = 1
                        items = line.strip().split("\t")[1:]
                        for index, value in enumerate(items):
                            if float(value) > 0:
                                dict_known[samples[index]] += 1
        with open(self.option("novel_mirna_count"), "r") as f1:
            head = f1.readline()
            samples = head.strip().split("\t")[1:]
            for sample in samples:
                dict_novel[sample] = 0
            for line in f1:
                items = line.strip().split("\t")[1:]
                for index, value in enumerate(items):
                    if float(value) > 0:
                        dict_novel[samples[index]] += 1
        with open(os.path.join(self.output_dir, "mirna_stat.xls"), "w") as w1:
            if self.option("known_mirna_count"):
                head = ["sample", "known_miRNAs", "novel_miRNAs", "total"]
                w1.write("\t".join(head) + "\n")
                for sample in self.samples:
                    w1.write(sample + "\t")
                    if sample in dict_known:
                        known = dict_known[sample]
                        w1.write(str(known) + "\t")
                    else:
                        known = 0
                        w1.write(str(known) + "\t")
                    if sample in dict_novel:
                        novel = dict_novel[sample]
                        w1.write(str(novel) + "\t")
                    else:
                        novel = 0
                        w1.write(str(novel) + "\t")
                    w1.write(str(int(known) + int(novel)) + "\n")
            else:
                head = ["sample", "number"]
                w1.write("\t".join(head) + "\n")
                for sample in self.samples:
                    w1.write(sample + "\t")
                    if sample in dict_novel:
                        novel = dict_novel[sample]
                        w1.write(str(novel) + "\n")
                    else:
                        novel = 0
                        w1.write(str(novel) + "\n")

        ## ncrna stat
        dict = {}
        self.logger.info("开始统计rfam结果")
        with open(self.option("rfam_summary"), "r") as f2:
            head = f2.readline()
            for line in f2:
                items = line.strip().split("\t")
                if items[0] in self.samples:
                    if items[0] not in dict:
                        dict[items[0]] = {}
                        dict[items[0]]["rRNA"] = items[1]
                        dict[items[0]]["tRNA"] = items[2]
                        dict[items[0]]["snoRNA"] = items[3]
                        dict[items[0]]["snRNA"] = items[4]
                    else:
                        dict[items[0]]["rRNA"] = items[1]
                        dict[items[0]]["tRNA"] = items[2]
                        dict[items[0]]["snoRNA"] = items[3]
                        dict[items[0]]["snRNA"] = items[4]
        self.logger.info("开始统计repeat结果")
        with open(self.option("repeat_summary"), "r") as f3:
            head = f3.readline()
            for line in f3:
                items = line.strip().split("\t")
                if items[0] in self.samples:
                    if items[0] not in dict:
                        dict[items[0]] = {}
                        dict[items[0]]["repbase"] = sum(map(int, items[1:]))
                    else:
                        dict[items[0]]["repbase"] = sum(map(int, items[1:]))
        self.logger.info("开始统计exon结果")
        with open(self.option("exon_summary"), "r") as f4:
            head = f4.readline()
            for line in f4:
                items = line.strip().split("\t")
                if items[0] in self.samples:
                    if items[0] not in dict:
                        dict[items[0]] = {}
                        dict[items[0]]["exon"] = items[1]
                    else:
                        dict[items[0]]["exon"] = items[1]
        self.logger.info("开始统计intron结果")
        with open(self.option("intron_summary"), "r") as f5:
            head = f5.readline()
            for line in f5:
                items = line.strip().split("\t")
                if items[0] in self.samples:
                    if items[0] not in dict:
                        dict[items[0]] = {}
                        dict[items[0]]["intron"] = items[1]
                    else:
                        dict[items[0]]["intron"] = items[1]
        self.logger.info("开始统计已知miRNA结果")
        if self.option("known_mrd"):
            with open(self.option("known_mrd"), "r") as f1:
                dict_mrd = {}
                for line in f1:
                    if re.compile(r'\S+_\d+_x\d+\s+').match(line):
                        if not line.split(" ")[0] in dict_mrd:
                            dict_mrd[line.split(" ")[0]] = 1
                            tmp_name = line.split(" ")[0].split("_")[0]
                            sample = self.parse_config(file=self.option("config").prop["path"], section="NAME", name=tmp_name)
                            number = int(line.split(" ")[0].split("_x")[1])
                            if sample not in dict:
                                dict[sample] = {}
                                dict[sample]["known_mirna"] = number
                            else:
                                if "known_mirna" not in dict[sample]:
                                    dict[sample]["known_mirna"] = number
                                else:
                                    dict[sample]["known_mirna"] += number
        self.logger.info("开始统计新预测miRNA结果")
        if self.option("method").lower() == "mirdeep2":
            with open(self.option("novel_mrd"), "r") as f1:
                dict_mrd = {}
                for line in f1:
                    if re.compile(r'\S+_\d+_x\d+\s+').match(line):
                        if not line.split(" ")[0] in dict_mrd:
                            dict_mrd[line.split(" ")[0]] = 1
                            tmp_name = line.split(" ")[0].split("_")[0]
                            sample = self.parse_config(file=self.option("config").prop["path"], section="NAME", name=tmp_name)
                            number = int(line.split(" ")[0].split("_x")[1])
                            if sample not in dict:
                                dict[sample] = {}
                                dict[sample]["novel_mirna"] = number
                            else:
                                if "novel_mirna" not in dict[sample]:
                                    dict[sample]["novel_mirna"] = number
                                else:
                                    dict[sample]["novel_mirna"] += number
        elif self.option("method").lower() == "mireap":
            with open(self.option("novel_mrd"), "r") as f1:
                dict_mrd = {}
                for line in f1:
                    if len(line.strip().split(" ")) >= 3:
                        mirna_id = line.strip().split(" ")[1]
                        if re.compile(r'^S.*_\d+$').match(mirna_id):
                            number = int(line.strip().split(" ")[2])
                            if mirna_id not in dict_mrd:
                                dict_mrd[mirna_id] = 1
                                tmp_name = mirna_id.split("_")[0]
                                sample = self.parse_config(file=self.option("config").prop["path"], section="NAME", name=tmp_name)
                                if sample not in dict:
                                    dict[sample] = {}
                                    dict[sample]["novel_mirna"] = number
                                else:
                                    if "novel_mirna" not in dict[sample]:
                                        dict[sample]["novel_mirna"] = number
                                    else:
                                        dict[sample]["novel_mirna"] += number
        else:
            novel_mirna_count_pd = pd.read_table(self.option("novel_mirna_count"), index_col=0)
            for sample in novel_mirna_count_pd.columns:
                sum_sample = novel_mirna_count_pd[sample].sum()
                if len(novel_mirna_count_pd[sample]) == 0:
                    sum_sample = 0
                if sample not in dict:
                    dict[sample] = {}
                dict[sample]["novel_mirna"] = sum_sample
        self.logger.info("开始统计所有miRNA序列结果")
        if self.option("input_type") == "raw":
            qc_stat_dir = os.path.join(self.option("qc_output"), "clean_data")
            for sample in self.samples:
                sample_stat = qc_stat_dir + "/" + sample + "_qc_stat.xls"
                if os.path.exists(sample_stat):
                    with open(sample_stat, "r") as f:
                        head = f.readline()
                        info = f.readline()
                        if sample not in dict:
                            dict[sample] = {}
                            dict[sample]["total"] = info.strip().split("\t")[6]
                        else:
                            dict[sample]["total"] = info.strip().split("\t")[6]
                else:
                    if sample not in dict:
                        dict[sample] = {}
                        dict[sample]["total"] = 0
                    else:
                        dict[sample]["total"] = 0
                print sample + "\t" + str(dict[sample]["total"]) + "\n"
        else:
            config = self.parse_config2(file=self.option("config").prop["path"])
            for seq_record in SeqIO.parse(self.option('input_fa').prop["path"], "fasta"):
                tmp_name = seq_record.id.split("_")[0]
                sample = config.get("NAME", tmp_name)
                number = int(seq_record.id.split("_x")[1])
                if sample not in dict:
                    dict[sample] = {}
                    dict[sample]["total"] = number
                else:
                    if "total" not in dict[sample]:
                        dict[sample]["total"] = number
                    else:
                        dict[sample]["total"] += number
        self.logger.info("开始写ncrna_stat.xls")
        with open(os.path.join(self.output_dir, "ncrna_stat.xls"), "w") as w1:
            head = ["Types"]
            head.extend(self.samples)
            w1.write("\t".join(head) + "\n")
            w1.write("rRNA" + "\t")
            for sample in self.samples[0:-1]:
                if sample in dict:
                    if "rRNA" in dict[sample] and int(dict[sample]["total"]) != 0:
                        w1.write(str(dict[sample]["rRNA"]) + "(" + str(round(float(dict[sample]["rRNA"])/float(dict[sample]["total"]), 4) * 100) + "%)" + "\t")
                    else:
                        w1.write("0(0.0%)" + "\t")
                else:
                    w1.write("0(0.0%)" + "\t")
            if self.samples[-1] in dict:
                if "rRNA" in dict[self.samples[-1]] and int(dict[self.samples[-1]]["total"]) != 0:
                    w1.write(str(dict[self.samples[-1]]["rRNA"]) + "(" + str(round(float(dict[self.samples[-1]]["rRNA"])/float(dict[self.samples[-1]]["total"]), 4) * 100) + "%)")
                else:
                        w1.write("0(0.0%)")
            else:
                w1.write("0(0.0%)")
            w1.write("\n")
            w1.write("tRNA" + "\t")
            for sample in self.samples[0:-1]:
                if sample in dict:
                    if "tRNA" in dict[sample] and int(dict[sample]["total"]) != 0:
                        w1.write(str(dict[sample]["tRNA"]) + "(" + str(round(float(dict[sample]["tRNA"])/float(dict[sample]["total"]), 4) * 100) + "%)" + "\t")
                    else:
                        dict[sample]["tRNA"] = 0
                        w1.write("0(0.0%)" + "\t")
                else:
                    dict[sample]["tRNA"] = 0
                    w1.write("0(0.0%)" + "\t")
            if self.samples[-1] in dict:
                if "tRNA" in dict[self.samples[-1]] and int(dict[self.samples[-1]]["total"]) != 0:
                    w1.write(str(dict[self.samples[-1]]["tRNA"]) + "(" + str(round(float(dict[self.samples[-1]]["tRNA"])/float(dict[self.samples[-1]]["total"]), 4) * 100) + "%)")
                else:
                    dict[self.samples[-1]]["tRNA"] = 0
                    w1.write("0(0.0%)")
            else:
                dict[self.samples[-1]]["tRNA"] = 0
                w1.write("0(0.0%)")
            w1.write("\n")
            w1.write("snoRNA" + "\t")
            for sample in self.samples[0:-1]:
                if sample in dict:
                    if "snoRNA" in dict[sample] and int(dict[sample]["total"]) != 0:
                        w1.write(str(dict[sample]["snoRNA"]) + "(" + str(round(float(dict[sample]["snoRNA"])/float(dict[sample]["total"]), 4) * 100) + "%)" + "\t")
                    else:
                        dict[sample]["snoRNA"] = 0
                        w1.write("0(0.0%)" + "\t")
                else:
                    dict[sample]["snoRNA"] = 0
                    w1.write("0(0.0%)" + "\t")
            if self.samples[-1] in dict:
                if "snoRNA" in dict[self.samples[-1]] and int(dict[self.samples[-1]]["total"]) != 0:
                    w1.write(str(dict[self.samples[-1]]["snoRNA"]) + "(" + str(round(float(dict[self.samples[-1]]["snoRNA"])/float(dict[self.samples[-1]]["total"]), 4) * 100) + "%)")
                else:
                    dict[self.samples[-1]]["snoRNA"] = 0
                    w1.write("0(0.0%)")
            else:
                dict[self.samples[-1]]["snoRNA"] = 0
                w1.write("0(0.0%)")
            w1.write("\n")
            w1.write("snRNA" + "\t")
            for sample in self.samples[0:-1]:
                if sample in dict:
                    if "snRNA" in dict[sample] and int(dict[sample]["total"]) != 0:
                        w1.write(str(dict[sample]["snRNA"]) + "(" + str(round(float(dict[sample]["snRNA"])/float(dict[sample]["total"]), 4) * 100) + "%)" + "\t")
                    else:
                        dict[sample]["snRNA"] = 0
                        w1.write("0(0.0%)" + "\t")
                else:
                    dict[sample]["snRNA"] = 0
                    w1.write("0(0.0%)" + "\t")
            if self.samples[-1] in dict:
                if "snRNA" in dict[self.samples[-1]] and int(dict[self.samples[-1]]["total"]) != 0:
                    w1.write(str(dict[self.samples[-1]]["snRNA"]) + "(" + str(round(float(dict[self.samples[-1]]["snRNA"])/float(dict[self.samples[-1]]["total"]), 4) * 100) + "%)")
                else:
                    dict[self.samples[-1]]["snRNA"] = 0
                    w1.write("0(0.0%)")
            else:
                dict[self.samples[-1]]["snRNA"] = 0
                w1.write("0(0.0%)")
            w1.write("\n")
            w1.write("repbase" + "\t")
            for sample in self.samples[0:-1]:
                if sample in dict:
                    if "repbase" in dict[sample] and int(dict[sample]["total"]) != 0:
                        w1.write(str(dict[sample]["repbase"]) + "(" + str(round(float(dict[sample]["repbase"])/float(dict[sample]["total"]), 4) * 100) + "%)" + "\t")
                    else:
                        dict[sample]["repbase"] = 0
                        w1.write("0(0.0%)" + "\t")
                else:
                    dict[sample]["repbase"] = 0
                    w1.write("0(0.0%)" + "\t")
            if self.samples[-1] in dict:
                if "repbase" in dict[self.samples[-1]] and int(dict[self.samples[-1]]["total"]) != 0:
                    w1.write(str(dict[self.samples[-1]]["repbase"]) + "(" + str(round(float(dict[self.samples[-1]]["repbase"])/float(dict[self.samples[-1]]["total"]), 4) * 100) + "%)")
                else:
                    dict[self.samples[-1]]["repbase"] = 0
                    w1.write("0(0.0%)")
            else:
                dict[self.samples[-1]]["repbase"] = 0
                w1.write("0(0.0%)")
            w1.write("\n")
            w1.write("exon" + "\t")
            for sample in self.samples[0:-1]:
                if sample in dict:
                    if "exon" in dict[sample] and int(dict[sample]["total"]) != 0:
                        w1.write(str(dict[sample]["exon"]) + "(" + str(round(float(dict[sample]["exon"])/float(dict[sample]["total"]), 4) * 100) + "%)" + "\t")
                    else:
                        dict[sample]["exon"] = 0
                        w1.write("0(0.0%)" + "\t")
                else:
                    dict[sample]["exon"] = 0
                    w1.write("0(0.0%)" + "\t")
            if self.samples[-1] in dict:
                if "exon" in dict[self.samples[-1]] and int(dict[self.samples[-1]]["total"]) != 0:
                    w1.write(str(dict[self.samples[-1]]["exon"]) + "(" + str(round(float(dict[self.samples[-1]]["exon"])/float(dict[self.samples[-1]]["total"]), 4) * 100) + "%)")
                else:
                    dict[self.samples[-1]]["exon"] = 0
                    w1.write("0(0.0%)")
            else:
                dict[self.samples[-1]]["exon"] = 0
                w1.write("0(0.0%)")
            w1.write("\n")
            w1.write("intron" + "\t")
            for sample in self.samples[0:-1]:
                if sample in dict:
                    if "intron" in dict[sample] and int(dict[sample]["total"]) != 0:
                        w1.write(str(dict[sample]["intron"]) + "(" + str(round(float(dict[sample]["intron"])/float(dict[sample]["total"]), 4) * 100) + "%)" + "\t")
                    else:
                        dict[sample]["intron"] = 0
                        w1.write("0(0.0%)" + "\t")
                else:
                    dict[sample]["intron"] = 0
                    w1.write("0(0.0%)" + "\t")
            if self.samples[-1] in dict:
                if "intron" in dict[self.samples[-1]] and int(dict[self.samples[-1]]["total"]) != 0:
                    w1.write(str(dict[self.samples[-1]]["intron"]) + "(" + str(round(float(dict[self.samples[-1]]["intron"])/float(dict[self.samples[-1]]["total"]), 4) * 100) + "%)")
                else:
                    dict[self.samples[-1]]["intron"] = 0
                    w1.write("0(0.0%)")
            else:
                dict[self.samples[-1]]["intron"] = 0
                w1.write("0(0.0%)")
            w1.write("\n")
            w1.write("unannotated" + "\t")
            for sample in self.samples[0:-1]:
                if sample in dict and int(dict[sample]["total"]) != 0:
                    if "known_mirna" not in dict[sample]:
                        dict[sample]["known_mirna"] = 0
                    dict[sample]["unannotated"] = int(dict[sample]["total"]) - int(dict[sample]["rRNA"]) - int(dict[sample]["tRNA"]) \
                                                  - int(dict[sample]["snoRNA"]) - int(dict[sample]["snRNA"]) - int(dict[sample]["repbase"]) \
                                                  - int(dict[sample]["exon"]) - int(dict[sample]["intron"]) - int(dict[sample]["known_mirna"])
                    w1.write(str(dict[sample]["unannotated"]) + "(" + str(round(float(dict[sample]["unannotated"])/float(dict[sample]["total"]), 4) * 100) + "%)" + "\t")
                else:
                    dict[sample]["unannotated"] = 0
                    w1.write("0(0.0%)" + "\t")
            if self.samples[-1] in dict and int(dict[self.samples[-1]]["total"]) != 0:
                if "known_mirna" not in dict[samples[-1]]:
                        dict[samples[-1]]["known_mirna"] = 0
                dict[self.samples[-1]]["unannotated"] = int(dict[samples[-1]]["total"]) - int(dict[samples[-1]]["rRNA"]) - int(dict[samples[-1]]["tRNA"]) \
                                                  - int(dict[samples[-1]]["snoRNA"]) - int(dict[samples[-1]]["snRNA"]) - int(dict[samples[-1]]["repbase"]) \
                                                  - int(dict[samples[-1]]["exon"]) - int(dict[samples[-1]]["intron"]) - int(dict[samples[-1]]["known_mirna"])
                w1.write(str(dict[self.samples[-1]]["unannotated"]) + "(" + str(round(float(dict[self.samples[-1]]["unannotated"])/float(dict[self.samples[-1]]["total"]), 4) * 100) + "%)")
            else:
                dict[samples[-1]]["unannotated"] = 0
                w1.write("0(0.0%)")
            w1.write("\n")
            w1.write("total" + "\t")
            for sample in self.samples[0:-1]:
                if sample in dict:
                    if "total" in dict[sample]:
                        w1.write(str(dict[sample]["total"]) + "\t")
                    else:
                        w1.write("0")
                else:
                    w1.write("0")
            if self.samples[-1] in dict:
                if "total" in dict[self.samples[-1]]:
                    w1.write(str(dict[self.samples[-1]]["total"]))
                else:
                        w1.write("0(0.0%)")
            else:
                w1.write("0(0.0%)")
            w1.write("\n")
        self.logger.info("开始写srna_stat.xls")
        with open(os.path.join(self.output_dir, "srna_stat.xls"), "w") as w1:
            head = ["Types"]
            head.extend(self.samples)
            w1.write("\t".join(head) + "\n")
            if self.option("known_mirna_count"):
                w1.write("known_mirna" + "\t")
                for sample in self.samples[0:-1]:
                    if sample in dict:
                        if "known_mirna" in dict[sample] and int(dict[sample]["total"]) != 0:
                            w1.write(str(dict[sample]["known_mirna"]) + "(" + str(round(float(dict[sample]["known_mirna"])/float(dict[sample]["total"]), 4) * 100) + "%)" + "\t")
                        else:
                            dict[sample]["known_mirna"] = 0
                            w1.write("0(0.0%)" + "\t")
                    else:
                        w1.write("0(0.0%)" + "\t")
                if self.samples[-1] in dict:
                    if "known_mirna" in dict[self.samples[-1]] and int(dict[self.samples[-1]]["total"]) != 0:
                        w1.write(str(dict[self.samples[-1]]["known_mirna"]) + "(" + str(round(float(dict[self.samples[-1]]["known_mirna"])/float(dict[self.samples[-1]]["total"]), 4) * 100) + "%)")
                    else:
                        dict[self.samples[-1]]["known_mirna"] = 0
                        w1.write("0(0.0%)")
                else:
                    w1.write("0(0.0%)")
                w1.write("\n")
            else:
                for sample in self.samples[0:-1]:
                    dict[sample]["known_mirna"] = 0
            w1.write("novel_mirna" + "\t")
            for sample in self.samples[0:-1]:
                if sample in dict:
                    if "novel_mirna" in dict[sample] and int(dict[sample]["total"]) != 0:
                        w1.write(str(dict[sample]["novel_mirna"]) + "(" + str(round(float(dict[sample]["novel_mirna"])/float(dict[sample]["total"]), 4) * 100) + "%)" + "\t")
                    else:
                        dict[sample]["novel_mirna"] = 0
                        w1.write("0(0.0%)" + "\t")
                else:
                    dict[sample]["novel_mirna"] = 0
                    w1.write("0(0.0%)" + "\t")
            if self.samples[-1] in dict:
                if "novel_mirna" in dict[self.samples[-1]] and int(dict[self.samples[-1]]["total"]) != 0:
                    w1.write(str(dict[self.samples[-1]]["novel_mirna"]) + "(" + str(round(float(dict[self.samples[-1]]["novel_mirna"])/float(dict[self.samples[-1]]["total"]), 4) * 100) + "%)")
                else:
                    dict[self.samples[-1]]["novel_mirna"] = 0
                    w1.write("0(0.0%)")
            else:
                dict[self.samples[-1]]["novel_mirna"] = 0
                w1.write("0(0.0%)")
            w1.write("\n")
            w1.write("rRNA" + "\t")
            for sample in self.samples[0:-1]:
                if sample in dict:
                    if "rRNA" in dict[sample] and int(dict[sample]["total"]) != 0:
                        w1.write(str(dict[sample]["rRNA"]) + "(" + str(round(float(dict[sample]["rRNA"])/float(dict[sample]["total"]), 4) * 100) + "%)" + "\t")
                    else:
                        dict[sample]["rRNA"] = 0
                        w1.write("0(0.0%)" + "\t")
                else:
                    dict[sample]["rRNA"] = 0
                    w1.write("0(0.0%)" + "\t")
            if self.samples[-1] in dict:
                if "rRNA" in dict[self.samples[-1]] and int(dict[self.samples[-1]]["total"]) != 0:
                    w1.write(str(dict[self.samples[-1]]["rRNA"]) + "(" + str(round(float(dict[self.samples[-1]]["rRNA"])/float(dict[self.samples[-1]]["total"]), 4) * 100) + "%)")
                else:
                    dict[self.samples[-1]]["rRNA"] = 0
                    w1.write("0(0.0%)")
            else:
                dict[self.samples[-1]]["rRNA"] = 0
                w1.write("0(0.0%)")
            w1.write("\n")
            w1.write("tRNA" + "\t")
            for sample in self.samples[0:-1]:
                if sample in dict:
                    if "tRNA" in dict[sample] and int(dict[sample]["total"]) != 0:
                        w1.write(str(dict[sample]["tRNA"]) + "(" + str(round(float(dict[sample]["tRNA"])/float(dict[sample]["total"]), 4) * 100) + "%)" + "\t")
                    else:
                        dict[sample]["tRNA"] = 0
                        w1.write("0(0.0%)" + "\t")
                else:
                    dict[sample]["tRNA"] = 0
                    w1.write("0(0.0%)" + "\t")
            if self.samples[-1] in dict:
                if "tRNA" in dict[self.samples[-1]] and int(dict[self.samples[-1]]["total"]) != 0:
                    w1.write(str(dict[self.samples[-1]]["tRNA"]) + "(" + str(round(float(dict[self.samples[-1]]["tRNA"])/float(dict[self.samples[-1]]["total"]), 4) * 100) + "%)")
                else:
                    dict[self.samples[-1]]["tRNA"] = 0
                    w1.write("0(0.0%)")
            else:
                dict[self.samples[-1]]["tRNA"] = 0
                w1.write("0(0.0%)")
            w1.write("\n")
            w1.write("snoRNA" + "\t")
            for sample in self.samples[0:-1]:
                if sample in dict:
                    if "snoRNA" in dict[sample] and int(dict[sample]["total"]) != 0:
                        w1.write(str(dict[sample]["snoRNA"]) + "(" + str(round(float(dict[sample]["snoRNA"])/float(dict[sample]["total"]), 4) * 100) + "%)" + "\t")
                    else:
                        dict[sample]["snoRNA"] = 0
                        w1.write("0(0.0%)" + "\t")
                else:
                    dict[sample]["snoRNA"] = 0
                    w1.write("0(0.0%)" + "\t")
            if self.samples[-1] in dict:
                if "snoRNA" in dict[self.samples[-1]] and int(dict[self.samples[-1]]["total"]) != 0:
                    w1.write(str(dict[self.samples[-1]]["snoRNA"]) + "(" + str(round(float(dict[self.samples[-1]]["snoRNA"])/float(dict[self.samples[-1]]["total"]), 4) * 100) + "%)")
                else:
                    dict[self.samples[-1]]["snoRNA"] = 0
                    w1.write("0(0.0%)")
            else:
                dict[self.samples[-1]]["snoRNA"] = 0
                w1.write("0(0.0%)")
            w1.write("\n")
            w1.write("snRNA" + "\t")
            for sample in self.samples[0:-1]:
                if sample in dict:
                    if "snRNA" in dict[sample] and int(dict[sample]["total"]) != 0:
                        w1.write(str(dict[sample]["snRNA"]) + "(" + str(round(float(dict[sample]["snRNA"])/float(dict[sample]["total"]), 4) * 100) + "%)" + "\t")
                    else:
                        dict[sample]["snRNA"] = 0
                        w1.write("0(0.0%)" + "\t")
                else:
                    dict[sample]["snRNA"] = 0
                    w1.write("0(0.0%)" + "\t")
            if self.samples[-1] in dict:
                if "snRNA" in dict[self.samples[-1]] and int(dict[self.samples[-1]]["total"]) != 0:
                    w1.write(str(dict[self.samples[-1]]["snRNA"]) + "(" + str(round(float(dict[self.samples[-1]]["snRNA"])/float(dict[self.samples[-1]]["total"]), 4) * 100) + "%)")
                else:
                    dict[self.samples[-1]]["snRNA"] = 0
                    w1.write("0(0.0%)")
            else:
                dict[self.samples[-1]]["snRNA"] = 0
                w1.write("0(0.0%)")
            w1.write("\n")
            w1.write("repbase" + "\t")
            for sample in self.samples[0:-1]:
                if sample in dict:
                    if "repbase" in dict[sample] and int(dict[sample]["total"]) != 0:
                        w1.write(str(dict[sample]["repbase"]) + "(" + str(round(float(dict[sample]["repbase"])/float(dict[sample]["total"]), 4) * 100) + "%)" + "\t")
                    else:
                        dict[sample]["repbase"] = 0
                        w1.write("0(0.0%)" + "\t")
                else:
                    dict[sample]["repbase"] = 0
                    w1.write("0(0.0%)" + "\t")
            if self.samples[-1] in dict:
                if "repbase" in dict[self.samples[-1]] and int(dict[self.samples[-1]]["total"]) != 0:
                    w1.write(str(dict[self.samples[-1]]["repbase"]) + "(" + str(round(float(dict[self.samples[-1]]["repbase"])/float(dict[self.samples[-1]]["total"]), 4) * 100) + "%)")
                else:
                    dict[self.samples[-1]]["repbase"] = 0
                    w1.write("0(0.0%)")
            else:
                dict[self.samples[-1]]["repbase"] = 0
                w1.write("0(0.0%)")
            w1.write("\n")
            w1.write("exon" + "\t")
            for sample in self.samples[0:-1]:
                if sample in dict:
                    if "exon" in dict[sample] and int(dict[sample]["total"]) != 0:
                        w1.write(str(dict[sample]["exon"]) + "(" + str(round(float(dict[sample]["exon"])/float(dict[sample]["total"]), 4) * 100) + "%)" + "\t")
                    else:
                        dict[sample]["exon"] = 0
                        w1.write("0(0.0%)" + "\t")
                else:
                    dict[sample]["exon"] = 0
                    w1.write("0(0.0%)" + "\t")
            if self.samples[-1] in dict:
                if "exon" in dict[self.samples[-1]] and int(dict[self.samples[-1]]["total"]) != 0:
                    w1.write(str(dict[self.samples[-1]]["exon"]) + "(" + str(round(float(dict[self.samples[-1]]["exon"])/float(dict[self.samples[-1]]["total"]), 4) * 100) + "%)")
                else:
                    dict[self.samples[-1]]["exon"] = 0
                    w1.write("0(0.0%)")
            else:
                dict[self.samples[-1]]["exon"] = 0
                w1.write("0(0.0%)")
            w1.write("\n")
            w1.write("intron" + "\t")
            for sample in self.samples[0:-1]:
                if sample in dict:
                    if "intron" in dict[sample] and int(dict[sample]["total"]) != 0:
                        w1.write(str(dict[sample]["intron"]) + "(" + str(round(float(dict[sample]["intron"])/float(dict[sample]["total"]), 4) * 100) + "%)" + "\t")
                    else:
                        dict[sample]["intron"] = 0
                        w1.write("0(0.0%)" + "\t")
                else:
                    dict[sample]["intron"] = 0
                    w1.write("0(0.0%)" + "\t")
            if self.samples[-1] in dict:
                if "intron" in dict[self.samples[-1]] and int(dict[self.samples[-1]]["total"]) != 0:
                    w1.write(str(dict[self.samples[-1]]["intron"]) + "(" + str(round(float(dict[self.samples[-1]]["intron"])/float(dict[self.samples[-1]]["total"]), 4) * 100) + "%)")
                else:
                    dict[self.samples[-1]]["intron"] = 0
                    w1.write("0(0.0%)")
            else:
                dict[self.samples[-1]]["intron"] = 0
                w1.write("0(0.0%)")
            w1.write("\n")
            w1.write("unknown" + "\t")
            for sample in self.samples[0:-1]:
                if sample in dict and int(dict[sample]["total"]) != 0:
                    dict[sample]["unknown"] = int(dict[sample]["total"]) - int(dict[sample]["rRNA"]) - int(dict[sample]["tRNA"]) \
                                                  - int(dict[sample]["snoRNA"]) - int(dict[sample]["snRNA"]) - int(dict[sample]["repbase"]) \
                                                  - int(dict[sample]["exon"]) - int(dict[sample]["intron"]) - int(dict[sample]["known_mirna"]) - int(dict[sample]["novel_mirna"])
                    w1.write(str(dict[sample]["unknown"]) + "(" + str(round(float(dict[sample]["unknown"])/float(dict[sample]["total"]), 4) * 100) + "%)" + "\t")
                else:
                    dict[sample]["unknown"] = 0
                    w1.write("0(0.0%)" + "\t")
            if self.samples[-1] in dict and int(dict[self.samples[-1]]["total"]) != 0:
                dict[self.samples[-1]]["unknown"] = int(dict[self.samples[-1]]["total"]) - int(dict[self.samples[-1]]["rRNA"]) - int(dict[self.samples[-1]]["tRNA"]) \
                                                  - int(dict[self.samples[-1]]["snoRNA"]) - int(dict[self.samples[-1]]["snRNA"]) - int(dict[self.samples[-1]]["repbase"]) \
                                                  - int(dict[self.samples[-1]]["exon"]) - int(dict[self.samples[-1]]["intron"]) - int(dict[self.samples[-1]]["known_mirna"]) - int(dict[self.samples[-1]]["novel_mirna"])
                w1.write(str(dict[self.samples[-1]]["unknown"]) + "(" + str(round(float(dict[self.samples[-1]]["unknown"])/float(dict[self.samples[-1]]["total"]), 4) * 100) + "%)")
            else:
                dict[self.samples[-1]]["unknown"] = 0
                w1.write("0(0.0%)")
            w1.write("\n")
            w1.write("total" + "\t")
            for sample in self.samples[0:-1]:
                if sample in dict:
                    if "total" in dict[sample]:
                        w1.write(str(dict[sample]["total"]) + "\t")
                    else:
                        w1.write("0" + "\t")
                else:
                    w1.write("0" + "\t")
            if self.samples[-1] in dict:
                if "total" in dict[self.samples[-1]]:
                    w1.write(str(dict[self.samples[-1]]["total"]))
                else:
                    w1.write("0")
            else:
                w1.write("0")
            w1.write("\n")

        self.logger.info("开始写srna_stat_for_graph.xls")
        if self.option("known_mirna_count"):
            with open(os.path.join(self.output_dir, "srna_stat_for_graph.xls"), "w") as w1:
                head = ["sample", "known_mirna", "novel_mirna", "rRNA", "snoRNA", "snRNA", "tRNA", "repbase", "exon", "intron", "unknow"]
                w1.write("\t".join(head) + "\n")
                for sample in self.samples:
                    self.logger.info(dict[sample])
                    if sample in dict:
                        w1.write("\t".join([sample, str(dict[sample]["known_mirna"]), str(dict[sample]["novel_mirna"]), str(dict[sample]["rRNA"]),
                                           str(dict[sample]["snoRNA"]), str(dict[sample]["snRNA"]), str(dict[sample]["tRNA"]), str(dict[sample]["repbase"]),
                                           str(dict[sample]["exon"]), str(dict[sample]["intron"]), str(dict[sample]["unknown"])]) + "\n")
        else:
            with open(os.path.join(self.output_dir, "srna_stat_for_graph.xls"), "w") as w1:
                head = ["sample", "novel_mirna", "rRNA", "snoRNA", "snRNA", "tRNA", "repbase", "exon",
                        "intron", "unknow"]
                w1.write("\t".join(head) + "\n")
                for sample in self.samples:
                    self.logger.info(dict[sample])
                    if sample in dict:
                        w1.write("\t".join([sample, str(dict[sample]["novel_mirna"]),
                                            str(dict[sample]["rRNA"]),
                                            str(dict[sample]["snoRNA"]), str(dict[sample]["snRNA"]),
                                            str(dict[sample]["tRNA"]), str(dict[sample]["repbase"]),
                                            str(dict[sample]["exon"]), str(dict[sample]["intron"]),
                                            str(dict[sample]["unknown"])]) + "\n")

class TestFunction(unittest.TestCase):

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "SrnaStat" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "small_rna_v2.srna.srna_stat",
            "instant": False,
            "options": dict(
                #known_mirna_count="/mnt/ilustre/users/sanger-dev/workspace/20181107/ath_srna18-09-12/Srna/output/known_mirna/known_mirna_count.xls",
                novel_mirna_count="/mnt/ilustre/users/sanger-dev/workspace/20200721/Single_srna09-16-31/Srna/Mireap/output/novel_mirna_count.xls",
                rfam_summary="/mnt/ilustre/users/sanger-dev/workspace/20200721/Single_srna09-16-31/Srna/RfamStat/output/rfam_summary.xls",
                repeat_summary="/mnt/ilustre/users/sanger-dev/workspace/20200721/Single_srna09-16-31/Srna/RepeatStat/output/repeat_summary.xls",
                exon_summary="/mnt/ilustre/users/sanger-dev/workspace/20200721/Single_srna09-16-31/Srna/GenomeStat/output/exon_summary.xls",
                intron_summary="/mnt/ilustre/users/sanger-dev/workspace/20200721/Single_srna09-16-31/Srna/GenomeStat/output/intron_summary.xls",
                input_fa="/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/miRNA/sRNA_test/small.fa",
                list="/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/miRNA/sRNA_test/list.txt",
                config="/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/miRNA/sRNA_test/qc_file.config",
                novel_mrd="/mnt/ilustre/users/sanger-dev/workspace/20200721/Single_srna09-16-31/Srna/KnownMirna/output/known_mirna.mrd",
                method="mireap"
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    unittest.main()

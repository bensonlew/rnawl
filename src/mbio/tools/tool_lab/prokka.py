# -*- coding: utf-8 -*-
# __author__ = 'zzg'
# version 1.0
# last_modify: 2020.03.01

import os, time
import re, shutil
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from mbio.packages.bac_comp_genome.common_function import link_dir

class ProkkaAgent(Agent):
    """
    prokka小工具
    采用conda进行
    """
    def __init__(self, parent):
        super(ProkkaAgent, self).__init__(parent)
        options = [
            {"name": "fasta", "type": "infile", "format": "sequence.fasta"},# 序列文件
            {"name": "sample_name", "type": "string"},# 样本名
            {"name": "rfam", "type": "string","default": "false"},
            {"name": "species", "type": "string", "default": "Bacteria"},
            #{"name": "evalue", "type": "float", "default": 1e-5},
            #{"name": "max_target_seqs", "type": "float", "default": 10},
        ]
        self.add_option(options)
        #self._memory_increase_step = 20

    def check_options(self):
        if not self.option("fasta").is_set:
            raise OptionError("必须输入序列文件")

    def set_resource(self):
        self._cpu = 8
        self._memory = '40G'

    def end(self):
        super(ProkkaAgent, self).end()

class ProkkaTool(Tool):
    def __init__(self, config):
        super(ProkkaTool, self).__init__(config)
        self.conda = self.config.SOFTWARE_DIR + "/program/miniconda3/"
        self.set_environ(PATH = self.conda +"bin")

    def run_anno(self):
        """
        prokka注释
        :return:
        """
        if os.path.exists(self.work_dir + "/result"):
            os.removedirs(self.work_dir + "/result")
        cmd = "{} {} --compliant --outdir result --prefix {} --kingdom {} " .format("program/miniconda3/bin/prokka", self.option("fasta").prop['path'], self.option("sample_name"), self.option("species"), )
        if self.option("rfam") == "yes":
            cmd += " --rfam"
        self.logger.info(cmd)
        command = self.add_command("run_anno", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_anno运行完成！")
        else:
            self.set_error("run_anno运行运行出错!")

    def run_stat(self):
        """
        统计长度
        :return:
        """
        N_base = 0;gene=0;cds=0;misc=0;rrna=0;repeat=0;trna=0;tmrna=0;scaff_num=0
        with open(self.option("fasta").prop['path']) as f:
            for i in f.readlines():
                if i.startswith(">"):
                    pass
                else:
                    value = i.strip()
                    N_base += (value.count('n') + value.count('N'))
        with open(self.work_dir + "/result/" + self.option("sample_name") + ".gff") as v:
            all_data = v.readlines()
            for i in all_data:
                if i.startswith("##"):
                    if i.startswith("##sequence"):
                        scaff_num += 1
                elif i.startswith(">"):
                    break
                else:
                    if i.split("\t")[2] == "gene":
                        gene += (abs(int(i.split("\t")[4]) - int(i.split("\t")[3])) + 1)
                    elif i.split("\t")[2] == "CDS":
                        cds += (abs(int(i.split("\t")[4]) - int(i.split("\t")[3])) + 1)
                    elif i.split("\t")[2] == "misc_RNA":
                        misc += (abs(int(i.split("\t")[4]) - int(i.split("\t")[3])) + 1)
                    elif i.split("\t")[2] == "repeat_region":
                        repeat += (abs(int(i.split("\t")[4]) - int(i.split("\t")[3])) + 1)
                    elif i.split("\t")[2] == "rRNA":
                        rrna += (abs(int(i.split("\t")[4]) - int(i.split("\t")[3])) + 1)
                    elif i.split("\t")[2] == "tRNA":
                        trna += (abs(int(i.split("\t")[4]) - int(i.split("\t")[3])) + 1)
                    elif i.split("\t")[2] == "tmRNA":
                        tmrna += (abs(int(i.split("\t")[4]) - int(i.split("\t")[3])) + 1)
        with open(self.work_dir + "/result/" + self.option("sample_name") + ".txt") as g,open(self.work_dir + "/result/" + "stat","w") as t:
            all_data2 = g.readlines()
            for i in all_data2:
                if i.startswith("contigs"):
                    t.write(i)
                elif i.startswith("bases"):
                    t.write(i)
                    n_ratio = round( float(N_base)/float(i.split(":")[1].strip()) * 100,2)
                elif i.startswith("gene"):
                    t.write(i.strip() + "\t" + str(gene) + "\n")
                elif i.startswith("CDS"):
                    t.write(i.strip() + "\t" + str(cds) + "\n")
                elif i.startswith("misc_RNA"):
                    t.write(i.strip() + "\t" + str(misc) + "\n")
                elif i.startswith("rRNA"):
                    t.write(i.strip() + "\t" + str(rrna) + "\n")
                elif i.startswith("repeat_region"):
                    t.write(i.strip() + "\t" + str(repeat) + "\n")
                elif i.startswith("tRNA"):
                    t.write(i.strip() + "\t" + str(trna) + "\n")
                elif i.startswith("tmRNA"):
                    t.write(i.strip() + "\t" + str(tmrna) + "\n")
            t.write("N_base" + "\t" +  str(n_ratio) + "\n")
            t.write("sample_name" + "\t" + self.option("sample_name") + "\n")

    def set_output(self):
        """
        设置结果文件目录
        """
        if os.path.exists(self.output_dir + "/result"):
            os.removedirs(self.output_dir + "/result")
        link_dir(self.work_dir + "/result", self.output_dir + "/result")
        for file in os.listdir(self.output_dir + "/result"):
            if file.endswith(".err") or file.endswith(".log") or file.endswith(".tsv") or file.endswith(".fna") or file.endswith(".txt"):
                os.remove(self.output_dir + "/result/" + file)
        self.logger.info("生成结果文件完成")

    def run(self):
        """
        运行
        """
        super(ProkkaTool, self).run()
        self.run_anno()
        self.run_stat()
        self.set_output()
        self.end()
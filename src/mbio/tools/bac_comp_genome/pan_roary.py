# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'
# version 1.0
# last_modify: 20119.10.17

import os
import re, shutil
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from mbio.packages.bac_comp_genome.common_function import link_dir
import subprocess

class PanRoaryAgent(Agent):
    """
    细菌比较基因组采用Roary进行聚类
    采用conda进行
    """
    def __init__(self, parent):
        super(PanRoaryAgent, self).__init__(parent)
        options = [
            {"name": "infile_dir", "type": "infile", "format": "sequence.fasta_dir"},#所有样本的蛋白文件夹
            #{"name": "gff_dir", "type": "infile", "format": "bac_comp_genome.gff_dir"},#多有样本的gff文件夹
            #{"name": "genome_dir", "type": "infile", "format": "bac_comp_genome.fasta_dir"},#基因组的序列
            {"name": "identity", "type": "float", "default": 0.5},#Minimum alignment coverage for blastp
            {"name": "inflation", "type": "float", "default": 1.5},#mcl 膨胀因子范围（1-5）
            {"name": "coverage", "type": "float", "default": 0.5},  # 给出cdhit的参数coverage
            #{"name": "sample", "type": "string"} #所有参与分析的样本名称并且以"+"相连
        ]
        self.add_option(options)
        self._memory_increase_step = 50  # 每次重运行增加内存50G by qingchen.zhang

    def check_options(self):
        if not self.option("infile_dir").is_set:
            raise OptionError("必须输入所有蛋白文件的文件夹")
        # if self.option("gff_dir").is_set:
        #     raise OptionError("必须输入gff的文件夹!")
        # if not self.option("genome_dir").is_set:
        #     raise OptionError("必须输入基因组文件或完成的scaffold文件")

    def set_resource(self):
        self._cpu = 4
        infile_path = self.option("infile_dir").prop['path']
        number_list = os.listdir(infile_path)
        number = len(number_list) / 2 ##几十个样本就用几十G
        total_memory = number * 1 + 10
        self._memory = '{}G'.format(total_memory)

    def end(self):
        super(PanRoaryAgent, self).end()

class PanRoaryTool(Tool):
    def __init__(self, config):
        super(PanRoaryTool, self).__init__(config)
        self.activate = self.config.SOFTWARE_DIR + "/program/miniconda3/bin/activate"
        self.roary = self.config.SOFTWARE_DIR + "/bioinfo/compare_genome/software/roary" ##设置Roary的conda环境
        self.sh = self.config.PACKAGE_DIR + "/bac_comp_genome/roary.sh"
        self.shell_path = "/program/sh"
        self.faa = self.option("infile_dir").prop['path']
        # self.gff = self.option("gff_dir").prop["path"]
        # self.genome = self.option("genome_dir").prop["path"]
        self.get_gff3 = self.config.PACKAGE_DIR + "/bac_comp_genome/require_gff3_from_seq.py"
        self.deacivate = self.config.SOFTWARE_DIR + "/program/miniconda3/bin/deactivate"
        self.conda = self.config.SOFTWARE_DIR + "/program/miniconda3/etc/profile.d/conda.sh"

    def creat_shell(self, input_dir, identity, inflation,coverage):
        """
        修改因为路径转换导致shell不能加载的问题， 重写一个shell，将运行的命令直接写入文件里面去
        :return:
        """
        self.fix_sh = os.path.join(self.work_dir, "roary.sh")
        with open(self.fix_sh, "w") as w:
            w.write("#! /bin/bash\n")
            w.write("source {}\n".format(self.conda))
            w.write("source {} {}\n".format(self.activate, self.roary))
            w.write("roary -e --mafft -p 8 {} -z -i {} -iv {} -ic {}\n".format(input_dir, identity, inflation, coverage))
            w.write("conda deactivate\n")

    def run_roary(self):
        """
        运行Roary软件进行聚类计算
        :return:
        """
        input_dir = self.faa + "/*.faa"
        identity = int(self.option("identity") * 100)
        coverage = int(self.option("coverage") * 100)
        self.creat_shell(input_dir, identity, self.option("inflation"), coverage)
        cmd = "{} {}".format(self.shell_path, self.fix_sh)
        self.logger.info(cmd)
        command = self.add_command("run_roary", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_roary运行完成！")
        else:
            self.set_error("run_roary运行完成运行出错!")
        # try:
        #     subprocess.check_output(cmd, shell=True)
        #     self.logger.info("run_roary运行完成！")
        # except subprocess.CalledProcessError:
        #     self.set_error("run_roary运行完成运行出错!")

    def run_stat(self):
        """
        对roary运行出来的结果进行整理和统计
        """
        stat_infile = os.path.join(self.work_dir, "gene_presence_absence.Rtab")
        cluster_infile = os.path.join(self.work_dir, "gene_presence_absence.csv")
        cluster_outfile = os.path.join(self.work_dir, "homologues_cluster.xls")
        # cluster = {}
        all_samples_list = []
        ##根据stat表统计出每个cluster的样本个数和基因个数
        with open(stat_infile, 'r') as f:
            lines = f.readlines()
            header = lines[0].strip().split("\t")
            all_samples = header[1:]
            for sample in all_samples:
                new_sample = sample.strip(".faa")
                all_samples_list.append(new_sample)
            # for line in lines[1:]:
            #     line = line.strip().split("\t")
            #     cluster_name = line[0]
            #     for i in range(1, len(all_samples_list)+1):
            #         if int(line[i]) != 0:
            #             cluster[cluster_name]["gene"] += int(line[i])
            #             cluster[cluster_name]["sample"] += 1
        ##根据输入的cluster表生成最后cluster与样本名称对应关系的大表
        with open(cluster_infile, 'r') as f2, open(cluster_outfile, 'w') as w:
            all_samples_name = "\t".join(all_samples_list)
            w.write("Cluster_ID\tSample_number\tGene_number\t" + all_samples_name + "\n")
            cluster_num = 0
            lines = f2.readlines()
            for line in lines[1:]:
                cluster_num += 1
                line = line.strip().split(",")
                cluster_name = "CLUSTER" + str(cluster_num)
                sample_num = int(eval(line[3]))#因为得到的结果中含有双引号
                gene_num = int(eval(line[4]))
                total_sample_list = [] #结果表是固定的，前14列为固定的统计信息，15列开始为记录的样本的cluster信息
                for i in range(14, len(all_samples_list)+14):
                    sample = eval(line[i])
                    if sample !='':
                        new_sample = sample.split("\t")
                        sample_name = ','.join(new_sample)
                        sample_real = sample_name
                    else:
                        sample_real = "-"
                    total_sample_list.append(sample_real)
                total_all_sample = "\t".join(total_sample_list)
                w.write("{}\t{}\t{}\t{}\n".format(cluster_name, sample_num, gene_num, total_all_sample))

    def set_output(self):
        """
        设置结果文件目录
        """
        outfile = os.path.join(self.output_dir, "homologues_cluster.xls")
        if os.path.exists(outfile):
            os.remove(outfile)
        os.link(self.work_dir + "/homologues_cluster.xls", outfile)
        self.logger.info("生成结果文件完成")

    def run(self):
        """
        运行
        """
        super(PanRoaryTool, self).run()
        self.run_roary()
        self.run_stat()
        self.set_output()
        self.end()

# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
# version 1.0
# last_modify: 20119.05.09

import os
import re, shutil
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from mbio.packages.bac_comp_genome.common_function import link_dir

class PanPgapAgent(Agent):
    """
    细菌比较基因组PGAP
    """
    def __init__(self, parent):
        super(PanPgapAgent, self).__init__(parent)
        options = [
            #{"name": "strains", "type": "sting"},  # 物种以+连接
            {"name": "method", "type": "string", "default": "GF"},  # GF for GeneFamily method,  and MP for MultiParanoid method
            {"name": "infile_dir", "type": "infile", "format": "bac_comp_genome.input_dir"},  # 输入文件夹
            {"name": "score", "type": "int", "default": 40},#Minimum score in blast
            {"name": "evalue", "type": "string", "default": "1e-10"},#Maximal E-value in blastall
            {"name": "coverage", "type": "float", "default": 0.5},#Minimum alignment coverage for two homologous proteins
            {"name": "local", "type": "float", "default": 0.25},#Minimum local alignment overlap in MP method
            {"name": "global", "type": "float", "default": 0.5},#Minimum global alignment overlap in MP method
            {"name": "identity", "type": "float", "default": 0.5},#Minimum alignment indentity for two homologous proteins
            {"name": "inflation", "type": "float", "default": 1.5}, ###膨胀系数
        ]
        self.add_option(options)
        self._memory_increase_step = 50  # 每次重运行增加内存50G by qingchen.zhang

    def check_options(self):
        # if not self.option("strains"):
        #     raise OptionError("必须设置参数strains的!")
        if self.option("method") not in ['GF','MP']:
            raise OptionError("必须设置参数method的type：GF or MP!")
        if not self.option("infile_dir").is_set:
            raise OptionError("必须设置参数infile_dir文件夹!")

    def set_resource(self):
        self._cpu = 4
        infile_path = self.option("infile_dir").prop['path']
        number_list = os.listdir(infile_path)
        number = len(number_list) / 2 ##几十个样本就用几十G
        total_memory = number * 1 + 10
        self._memory = '{}G'.format(total_memory)

    def end(self):
        super(PanPgapAgent, self).end()

class PanPgapTool(Tool):
    def __init__(self, config):
        super(PanPgapTool, self).__init__(config)
        self.path = self.config.SOFTWARE_DIR + "/bioinfo/compare_genome/software/PGAP-1.2.1:" + self.config.SOFTWARE_DIR + "/bioinfo/compare_genome/software/mcl/bin:" + self.config.SOFTWARE_DIR + \
                    "/bioinfo/align/ncbi-blast-2.3.0+/bin:" + self.config.SOFTWARE_DIR + "bioinfo/phylogenetic/standard-RAxML-master:" + self.config.SOFTWARE_DIR + "/bioinfo/compare_genome/software/phylip-3.697/exe:"
        self.perl5path = self.config.SOFTWARE_DIR + "/bioinfo/compare_genome/software/PGAP-1.2.1"
        #self.strains = self.option("strains")
        self.fasta = self.option("infile_dir").prop['path']
        self.perl_path = "/program/perl/perls/perl-5.24.0/bin/perl"
        self.set_environ(PATH=self.path, PERL5LIB=self.perl5path)
        self.perl_script = self.config.SOFTWARE_DIR + "/bioinfo/compare_genome/software/PGAP-1.2.1/PGAP.pl"
        self.out = self.work_dir + '/out_result'
        self.format_convert = self.config.PACKAGE_DIR + "/tool_lab/format_convert.py"
        self.python = "/program/Python/bin/python"


    def format_covert(self):
        """
        对输入文件进行标准化处理，为序列名称加样本名称并以“|”分割，为变异分析整理成pgap的输入格式.pep .nuc .function
        .function文件为功能分析的文件，现转了一个只有标题的文件
        :return:
        """
        self.logger.info("开始进行格式转换")
        self.fasta_dir = self.work_dir + '/fasta_dir'
        if os.path.exists(self.fasta_dir):
            shutil.rmtree(self.fasta_dir)
        os.mkdir(self.fasta_dir)
        self.input_dir = self.work_dir + '/input_dir'
        if os.path.exists(self.input_dir):
            shutil.rmtree(self.input_dir)
        os.mkdir(self.input_dir)
        cmd = "{} {} -i {} -out {} -o {} -t 4".format(self.python, self.format_convert, self.fasta, self.fasta_dir, self.input_dir)
        self.logger.info(cmd)
        command = self.add_command("format_convert", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("格式转换运行完成！")
        else:
            self.set_error("格式转换运行出错!")

    def run_pagp(self):
        """
        运行pgap进行聚类
        :return:
        """
        self.get_strains()
        if os.path.exists(self.out):
            shutil.rmtree(self.out)
        os.mkdir(self.out)
        cmd = "{} {} --strains {} --cluster --input {} --output {} --cluster --identity {}  --thread 16 --coverage {}".format(self.perl_path, self.perl_script, self.strains, self.fasta_dir, self.out, self.option("identity"), self.option("coverage"))
        if self.option("method") in ["GF"]:
            cmd += " --method GF --score {} --evalue {} --inflation {}".format(self.option("score"), self.option("evalue"), self.option("inflation"))
        elif self.option("method") in ["MP"]:
            cmd += " --method MP --local {} --global {}".format(self.option("local"), self.option("global"))
        self.logger.info(cmd)
        command = self.add_command("run_pagp", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_pagp运行完成！")
        else:
            self.set_error("run_pagp运行完成运行出错!")

    def run_stat(self):
        """
        对PGAP的结果进行整理和统计
        :return:
        """
        self.logger.info("开始对聚类结果进行统计")
        cluster_result = os.path.join(self.out, "1.Orthologs_Cluster.txt")
        results = os.path.join(self.output_dir, "homologues_cluster.xls")
        with open(cluster_result, 'r') as f, open(results, 'w') as w:
            lines = f.readlines()
            header = lines[0]
            head_list = header.strip().split("\t")
            sample_list = head_list[1:]
            all_sample = '\t'.join(sample_list)
            w.write("Cluster_ID\tSample_number\tGene_number\t" + all_sample + "\n")
            cluster_num = 0
            all_sample_number = len(sample_list)+1
            for line in lines[1:]:
                cluster_num += 1
                line = line.strip().split("\t")
                sample_num = 0
                gene_num = 0
                for i in range(1, all_sample_number):
                    sample_name = line[i]
                    if sample_name != "-":
                        sample_num += 1
                    sample_cluster_list = sample_name.split(",")
                    for i in sample_cluster_list:
                        if i != "-":
                            gene_num += 1
                newline = "\t".join(line[1:])
                cluster_id = "CLUSTER" + str(cluster_num)
                w.write(cluster_id + "\t" + str(sample_num) + "\t" + str(gene_num) + "\t" + newline + "\n")
        self.logger.info("聚类结果统计完成")

    def get_strains(self):
        """
        根据传入的文件夹内的文件生成对应的样本名称
        :return:
        """
        all_files = os.listdir(self.fasta_dir)
        sample_list = []
        for file in all_files:
            if re.search(r".pep", file):
                sample = file.split(".pep")[0]
                sample_list.append(sample)
        self.strains = "+".join(sample_list)

    def set_output(self):
        """
        设置结果文件目录
        :return:
        """
        work_result = os.path.join(self.output_dir, 'homologues_cluster.xls')
        if not os.path.exists(work_result):
            self.set_error("不存在运行的结果目录文件")

    def run(self):
        super(PanPgapTool, self).run()
        self.format_covert()
        self.run_pagp()
        self.run_stat()
        self.set_output()
        self.end()

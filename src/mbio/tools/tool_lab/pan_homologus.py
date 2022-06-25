# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'
# version 1.0
# last_modify: 2019.1015

import os
import re, shutil
import pandas as pd
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from mbio.packages.bac_comp_genome.common_function import link_dir
from Bio import SeqIO

class PanHomologusAgent(Agent):
    """
    细菌比较基因组采用GET_HOMOLOGUS进行聚类分析，本质是blast+mcl的分析，主要是本软件还有很多操作
    """
    def __init__(self, parent):
        super(PanHomologusAgent, self).__init__(parent)
        options = [
            {"name": "method", "type": "string", "default": "omcl"},  # 三种方法BDBH、OMCL、OCOG
            {"name": "infile_dir", "type": "infile", "format": "sequence.fasta_dir"},  # 输入蛋白文件夹
            {"name": "coverage", "type": "float", "default": 0.5},#覆盖度blast的coverage
            {"name": "identity", "type": "float", "default": 0.5},#聚类的相似性值identity
            {"name": "evalue", "type": "string", "default": 1e-05},#聚类的相似性值
            {"name": "thread", "type": "int", "default": 4},  # 线程数和cpu数
            {"name": "filter", "type": "bool", "default": False}, #对聚类结果是否按长度进行过滤
            {"name": "filter_length", "type": "int", "default": 50}, #过滤长度
            {"name": "inflation", "type": "float", "default": 1.5}, #mcl膨胀因子
        ]
        self.add_option(options)
        self._memory_increase_step = 50  # 每次重运行增加内存50G by qingchen.zhang

    def check_options(self):
        """
        重写参数检查
        :return:
        """
        if not self.option("infile_dir").is_set:
            raise OptionError("必须设置参数infile_dir文件夹!")

    def set_resource(self):
        """
        设置资源使用情况
        :return:
        """
        self._cpu = 4
        infile_path = self.option("infile_dir").prop['path']
        number_list = os.listdir(infile_path)
        number = len(number_list) / 2 ##几十个样本就用几十G
        total_memory = number * 1 + 10
        self._memory = '{}G'.format(total_memory)

    def end(self):
        super(PanHomologusAgent, self).end()

class PanHomologusTool(Tool):
    def __init__(self, config):
        super(PanHomologusTool, self).__init__(config)
        self.path = self.config.SOFTWARE_DIR + "/bioinfo/compare_genome/software/get_homologues-master:" + self.config.SOFTWARE_DIR + "/program/R-3.3.1/bin:"
        self.perl5path = self.config.SOFTWARE_DIR + "/miniconda2/bin/perl"
        self.library = self.config.SOFTWARE_DIR + "/program/R-3.3.1/lib64/R/lib"
        self.r_home = self.config.SOFTWARE_DIR + "/program/R-3.3.1/lib64/R"
        self.c_include = self.config.SOFTWARE_DIR + "/program/R-3.3.1/lib64/R/include"
        self.set_environ(PATH=self.path, PERL5LIB=self.perl5path)
        self.set_environ(R_HOME=self.r_home, LD_LIBRARY_PATH=self.library, C_INCLUDE_PATH=self.c_include)
        self.fasta = self.option("infile_dir").prop['path']
        self.perl_path = "/miniconda2/bin/perl"
        self.perl_script = self.config.SOFTWARE_DIR + "/bioinfo/compare_genome/software/get_homologues-master/"
        self.out = self.work_dir + '/out_result'

    def check_and_get_faa(self):
        infile_path = self.option("infile_dir").prop['path']
        all_files = os.listdir(infile_path)
        for file in all_files:
            if file.endswith('.fna'):
                file_name = file.split('.fna')[0]
                fna_tmp = infile_path + '/' + file
                faa_tmp = infile_path+'/'+file_name+'.faa'
                if not os.path.exists(faa_tmp):
                    r_seq = SeqIO.parse(fna_tmp, "fasta")  #parse
                    with open(faa_tmp,'w') as fw:
                        for e in r_seq:
                            try:
                                faa_seq = e.seq.translate(table="Bacterial", cds=True)
                                fw.write('>'+e.description+'\n')
                                fw.write(str(faa_seq)+'\n')
                            except Exception as e:
                                print(e)

    def run_homologus(self,method):
        """
        运行GET_Homologus进行分析
        :return:
        """
        if self.option("coverage"):
            self.coverage = self.option("coverage")
        if self.option("identity"):
            self.identity = self.option("identity")
        self.logger.info("开始运行get_homologus")
        cmd = "{} {}get_homologues.pl -d {} -X -n {} -c -C {} -E {}".format(self.perl_path, self.perl_script, self.fasta, self.option("thread"), self.coverage, self.option("evalue"))
        if method in ["omcl"]:
            cmd += " -t 0 -M -S {} -F {}".format(self.identity, self.option("inflation"))
            cmd_name = "run_omcl"
        elif method in ["ocog"]:
            cmd += " -t 0 -G"
            cmd_name = "run_cog"
        else: #BDBH
            cmd += " -t 1 -S {}".format(self.identity)
            cmd_name = "run_bdbh"
        if self.option("filter"):
            cmd += " -f {}".format(self.option("filter_length"))
        self.logger.info(cmd)
        command = self.add_command(cmd_name, cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_homologus运行完成！")
        else:
            self.set_error("run_homologus运行完成运行出错!")

    def run_compare(self):
        """
        对第一步的结果进行挑选和不同方法之间的比较
        :return:
        """
        cmd = "{} {}compare_clusters.pl -o {} -m -T -d ".format(self.perl_path, self.perl_script, self.out)
        file_name = os.path.basename(self.fasta)
        real_name = file_name + "_homologues"
        dir_path = os.path.join(self.work_dir, real_name)
        self.logger.info(dir_path)
        method_list = []
        for files in os.listdir(dir_path):
            method_path = os.path.join(dir_path, files)
            if os.path.isdir(method_path):
                if re.search(r'OMCL', files):
                    method_list.append(method_path)
                elif re.search(r'COG', files):
                    method_list.append(method_path)
                elif re.search(r"BDBH", files):
                    method_list.append(method_path)
                # else:
                #     self.set_error("running error in homologus")
        method_path_name = ','.join(method_list)
        cmd += method_path_name
        self.logger.info(cmd)
        command = self.add_command("run_cluster", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_cluster运行完成！")
        else:
            self.set_error("run_cluster运行出错!")

    def get_cluster(self):
        """
        根据第二步产生的结果，对聚类结果进行处理
        :return:
        """
        file_name = os.path.basename(self.fasta)
        cluster_path = os.path.join(self.out, "pangenome_matrix_genes_t0.tab")
        data = pd.read_table(cluster_path, sep='\t', header=0)
        data1 = data.T
        data1 = data1.iloc[0:-1]
        data1 = data1.where(data1.notnull(), '-')
        data1.to_csv("all_sample.xls", sep='\t')
        os.system("sed -i '1d' all_sample.xls") #对第一行的表头删除转置的时候自动生成
        cluster_result = os.path.join(self.output_dir, "homologues_cluster.xls")
        with open("all_sample.xls", "r") as f, open(cluster_result, 'w') as w:
            header = f.readline()
            sample_list = []
            head_list = header.strip().split("\t")
            for sample in head_list[1:]:
                new_sample = sample.rstrip(".faa")
                sample_list.append(new_sample)# 获得所有的样本名称
            all_sample = "\t".join(sample_list)
            w.write("Cluster_ID\tSample_number\tGene_number\t" + all_sample + "\n")
            cluster_num = 0
            for line in f:
                cluster_num += 1
                total_sample_num = 0
                total_gene_num = 0
                line = line.strip().split("\t")
                total_sample_list = []
                sample_num = len(sample_list) + 1 #在样本数目上加1
                for i in range(1,sample_num):
                    new_line = line[i].rstrip(",")
                    if new_line != "-":
                        total_sample_num += 1
                    gene_name_list = new_line.split(",")
                    for j in gene_name_list:
                        if j != "-":
                            total_gene_num += 1
                    total_sample_list.append(new_line)
                cluster_name = "CLUSTER" + str(cluster_num)
                newline = "\t".join(total_sample_list)
                w.write(cluster_name +"\t" + str(total_sample_num) + "\t" + str(total_gene_num) + "\t" + newline + "\n")
        self.logger.info("完成统计和整理聚类结果")

    def set_output(self):
        """
        设置结果文件目录
        :return:
        """
        if os.path.exists(os.path.join(self.output_dir, "homologues_cluster.xls")):
            self.logger.info("生成了正确的结果文件")
        else:
            self.set_error("未能生成正确的结果文件")

    def run(self):
        super(PanHomologusTool, self).run()
        self.check_and_get_faa()
        methods = self.option("method").split(",")
        for method in methods:
            self.run_homologus(method)
        self.run_compare()
        self.get_cluster()
        self.set_output()
        self.end()

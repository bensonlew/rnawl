# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
# version 1.0
# last_modify: 20119.05.09

import os
import json
import re, shutil
import pandas as pd
from biocluster.agent import Agent
from biocluster.tool import Tool
from Bio import SeqIO
from biocluster.core.exceptions import OptionError
from mbio.packages.bac_comp_genome.common_function import link_dir

class PanOrthofinderAgent(Agent):
    """
    细菌比较基因组orthofinder软件进行计算
    """
    def __init__(self, parent):
        super(PanOrthofinderAgent, self).__init__(parent)
        options = [
            {"name": "thread", "type": "int", "default": 4},  # 线程数
            {"name": "infile_dir", "type": "infile", "format": "bac_comp_genome.input_dir"},  # 输入蛋白文件夹
            {"name": "program", "type": "string", "default": "diamond"},#Sequence search program [Default = diamond]")
            {"name": "inflation", "type": "float", "default": 1.5},#Minimum score in blast
            {"name": "species", "type": "bool", "default": False},#是否将物种名称加到序列id上
            {"name": "work_name", "type": "string", "default": "orthofinder"},#工作目录的名称
            {"name": "stop", "type": "string", "default": "og"},#工作流在工作完哪里之后就停下
            {"name": "identity", "type": "float", "default": 0.5},  ##给出cdhit的参数identity
            {"name": "coverage", "type": "float", "default": 0.5},  # 给出cdhit的参数coverage
        ]
        self.add_option(options)
        self._memory_increase_step = 50  # 每次重运行增加内存50G by qingchen.zhang

    def check_options(self):
        if not self.option("infile_dir").is_set:
            raise OptionError("必须设置参数infile_dir文件夹!")

    def set_resource(self):
        self._cpu = 4
        infile_path = self.option("infile_dir").prop['path']
        number_list = os.listdir(infile_path)
        number = len(number_list) / 2 ##几十个样本就用几十G
        total_memory = number * 1 + 10
        self._memory = '{}G'.format(str(total_memory))

    def end(self):
        super(PanOrthofinderAgent, self).end()

class PanOrthofinderTool(Tool):
    def __init__(self, config):
        super(PanOrthofinderTool, self).__init__(config)
        self.path = self.config.SOFTWARE_DIR + "/bioinfo/compare_genome/software/OrthoFinder-master/:" + self.config.SOFTWARE_DIR + "/bioinfo/compare_genome/software/mcl/bin:" + self.work_dir + "/orthofinder:"+ self.config.SOFTWARE_DIR + "/bioinfo/align/ncbi-blast-2.3.0+/bin:"+ self.config.SOFTWARE_DIR + "/bioinfo/align/diamond-0.8.35/"
        self.perl5path = self.config.SOFTWARE_DIR + "/program/perl/perls/perl-5.24.0/bin/perl"
        self.set_environ(PATH=self.path, PERL5LIB=self.perl5path)
        self.fasta = self.option("infile_dir").prop['path']
        self.python = "/program/Python/bin/python"
        self.python_script = self.config.SOFTWARE_DIR + "/bioinfo/compare_genome/software/OrthoFinder-master/orthofinder"
        self.out = self.work_dir + '/out_result'
        self.config_json = self.config.SOFTWARE_DIR + "/bioinfo/compare_genome/software/OrthoFinder-master/orthofinder/config.json"
        self.identity = int(self.option("identity")*100)
        self.coverage = int(self.option("coverage")*100)

    def get_config(self):
        """
        根据传入参数identity和coverage修改config文件
        :return:
        """
        self.orthofinder_path = os.path.join(self.work_dir, "orthofinder")
        if os.path.exists(self.orthofinder_path):
            shutil.rmtree(self.orthofinder_path)
        shutil.copytree(self.python_script, self.orthofinder_path)
        self.config_path = os.path.join(self.work_dir, 'orthofinder/config.json')
        if os.path.exists(self.config_path):
            os.remove(self.config_path)
        with open(self.config_json, "r") as f, open(self.config_path, 'w') as w:
            lines = f.readlines()
            for line in lines:
                if re.search(r'diamond blastp', line):
                    new_line_key = line.split(":")[0]
                    new_line_value = line.split(":")[1]
                    aa = eval(new_line_value)
                    bb = aa + " --id {} --query-cover {}".format(self.identity, self.coverage)
                    new_line = new_line_key + ":" + "\"{}\"".format(bb)
                    w.write(new_line + "\n")
                else:
                    w.write(line)

    def run_orthofinder(self):
        """
        运行orthofinder进行聚类或者拿到同源蛋白
        :return:
        """
        # self.get_config()
        if os.path.exists(self.out):
            shutil.rmtree(self.out)
        cmd = "{} {}/orthofinder.py -f {} -t {} -I {} -n {} -o {} -S {}".format(self.python, self.orthofinder_path, self.fasta,self.option("thread"), self.option("inflation"), self.option("work_name"), self.out, self.option("program"))
        if self.option("species"):
            cmd += " -X"
        elif self.option("stop"):
            cmd += " -{} ".format(self.option("stop"))
        self.logger.info(cmd)
        command = self.add_command("run_orthofinder", cmd).run()

        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_orthofinder运行完成！")
        else:
            self.set_error("run_orthofinder运行出错!")

    def make_table(self):
        """
        对orthofinder的结果进行整理
        :return:
        """
        self.logger.info("正在对orthofinder的结果进行整理")
        if os.path.exists(self.work_dir + '/Orthologs_Cluster1.xls'):
            os.remove(self.work_dir + '/Orthologs_Cluster1.xls')
        if os.path.exists(self.work_dir + '/Orthologs_Cluster.xls'):
            os.remove(self.work_dir + '/Orthologs_Cluster.xls')
        result_name = "Results_%s/Orthogroups" %(self.option("work_name"))
        result_path = os.path.join(self.out, result_name)
        cluster_result1 = os.path.join(result_path, 'Orthogroups.tsv') #paralogus
        cluster_result2 = os.path.join(result_path, 'Orthogroups_UnassignedGenes.tsv') #paralogus未比对上的结果
        os.system("head -n 1 %s > header.txt" %cluster_result1) ##保留表头
        os.system("cat %s >> Orthologs_Cluster.xls" %cluster_result1)
        os.system("cat %s >> Orthologs_Cluster.xls" %cluster_result2) ## 将结果1和2合并起来

        input_file = os.path.join(self.work_dir, "Orthologs_Cluster.xls")
        last_result = os.path.join(self.output_dir, 'homologues_cluster.xls')
        data = pd.read_table(input_file, sep='\t', header=0)
        data.fillna('-',inplace=True)
        data.to_csv('Orthologs_Cluster1.xls', index =0,sep='\t')
        new_input_file = os.path.join(self.work_dir, "Orthologs_Cluster1.xls")
        cluster_num = 0
        cluster = {}
        with open(last_result, 'w') as w,open(new_input_file, 'r') as f:
            lines = f.readlines()
            header = lines[0]
            all_samples = header.strip().split("\t")
            all_samples_name = "\t".join(all_samples[1:])
            w.write("Cluster_ID\tSample_number\tGene_number\t" + all_samples_name + "\n")
            for line in lines[1:]:
                line = line.strip().split("\t")
                total_sample_num = 0
                total_gene_num = 0
                if line[0] in ['Orthogroup']:
                    pass
                else:
                    cluster_num += 1
                    sample_num = len(all_samples)#在样本数目上加1
                    for i in range(1,sample_num):
                        # self.logger.info(line[0])
                        all_gene_name = line[i]
                        if all_gene_name != "-":
                            total_sample_num += 1
                        gene_name_list = line[i].split(",")
                        for i in gene_name_list:
                            if i != "-":
                                total_gene_num += 1
                    cluster[line[0]] = "CLUSTER" + str(cluster_num)
                    cluster_name = "CLUSTER" + str(cluster_num)
                    strip_space_line = [",".join([y.strip() for y in x.strip().split(",")]) for x in line[1:]]
                    newline = "\t".join(strip_space_line)
                    w.write(cluster_name +"\t" + str(total_sample_num) + "\t" + str(total_gene_num) + "\t" + newline + "\n")


    def get_fasta(self, cluster):
        """
        根据OG0000892与cluster对应关系，获得聚类得到的fasta结果
        :return:
        """
        fasta_name = "Results_%s/Orthogroup_Sequences" %(self.option("work_name"))
        fasta_path = os.path.join(self.out, fasta_name)
        cluster_fasta = os.path.join(self.output_dir, "/all_cluster.faa")
        with open(cluster_fasta, 'w') as w:
            for key in cluster.keys():
                file_name = key + '.fa'
                file_path = os.path.join(fasta_path, file_name)
                seq = ''
                for seq_record in SeqIO.parse(file_path, 'fasta'):#打开第一次主要是获取到聚类的结果id和样本名称
                    seq = seq_record.seq
                    break
                seq_id = cluster[key]
                w.write(">{}\n{}\n".format(seq_id,seq))

    def set_output(self):
        """
        设置结果目录
        :return:
        """
        self.logger.info("正在设置结果文件目录")
        if os.path.exists(os.path.join(self.output_dir, 'homologues_cluster.xls')):
            self.logger.info("生成了正确的聚类结果文件")
        if os.path.exists(os.path.join(self.output_dir, 'all_cluster.faa')):
            self.logger.info("生成了正确的聚类fasta序列文件")
        self.logger.info("生成结果文件完成")

    def run(self):
        """
        运行
        :return:
        """
        self.logger.info("开始运行tool！")
        super(PanOrthofinderTool, self).run()
        self.get_config()
        if not os.path.exists(self.out):
            self.run_orthofinder()
        self.make_table()
        #self.get_fasta(cluster)
        self.set_output()
        self.end()

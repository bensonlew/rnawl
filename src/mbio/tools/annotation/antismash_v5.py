# -*- coding: utf-8 -*-
# __author__ = 'zhaozhigang'
# last_modify: 2020.11.11
from biocluster.agent import Agent
from biocluster.tool import Tool
import os, re, shutil
from biocluster.core.exceptions import OptionError
import subprocess
import unittest
import pandas as pd


class AntismashV5Agent(Agent):
    """
    功能：antismash注释,即次级代谢分析
    输入：组装得到的基因组序列(必须是*.fasta后缀)或者gbk文件
    比对软件：antismash-5.1.2
    """

    def __init__(self, parent):
        super(AntismashV5Agent, self).__init__(parent)
        options = [
            {"name": "genome_seq", "type": "infile", "format": "sequence.fasta"},
            # 组装得到的基因组文件，尾缀必须是.fasta， 该步先进行基因预测（细菌用prodigal，真菌用glimmerhmm，否则报错）
            {"name": "genome_gbk", "type": "infile", "format": "gene_structure.gbk"},  # 样品基因组的gbk格式文件，无需进行预测 gene_structure.gbk
            {"name": "taxon", "type": "string", "default": "bacteria"},
            # 基因组序列来源：bacteria,fungi #目前软件中当输入为*.fasta，且taxon 为 fungi时，没有对应的基因预测方法
            {"name": "predict", "type": "string", "default": "none"},
            # 基因预测的方法，默认none，可选{glimmerhmm,prodigal,prodigal-m,none,error}
            {"name": "strictness", "type": "string", "default": "relaxed"}
            # 定义用于基于HMM的群集检测的严格程度（默认值：relaxed）strict严格，loose宽松
        ]
        self.add_option(options)
        #self.queue = 'BLAST'  # 投递到指定的队列BLAST

    def check_options(self):
        if not self.option("genome_seq").is_set and not self.option("genome_gbk").is_set:
            raise OptionError("必须设置参数genome_seq或者genome_gbk", code="31200101")
        return True

    def set_resource(self):
        self._cpu = 10
        self._memory = '70G'

    def end(self):
        super(AntismashV5Agent, self).end()


class AntismashV5Tool(Tool):
    def __init__(self, config):
        super(AntismashV5Tool, self).__init__(config)
        self.file_name = ""
        self.antismash = os.path.join(self.config.SOFTWARE_DIR, "bioinfo/Genomic/Sofware/antismash-5.1.2/")
        self.path = self.config.SOFTWARE_DIR + "/program/Python35/bin:" + self.config.SOFTWARE_DIR + '/program/sun_jdk1.8.0/bin:' + self.config.SOFTWARE_DIR + "/bioinfo/annotation/prodigal/Prodigal-2.6.3:" + self.config.SOFTWARE_DIR + "/bioinfo/align/diamond-0.8.36:" + self.config.SOFTWARE_DIR + "bioinfo/Genomic/Sofware/antismash-5.1.2:" + self.config.SOFTWARE_DIR + "/bioinfo/Genomic/Sofware/glimmer3.02/bin:" + self.config.SOFTWARE_DIR + "/bioinfo/align/meme/meme/bin/:" + self.config.SOFTWARE_DIR + "/bioinfo/align/hmmer-3.1b2-linux-intel-x86_64/binaries:" + self.config.SOFTWARE_DIR + "/bioinfo/Genomic/Sofware/hmmer-2.3.2/bin:" + self.config.SOFTWARE_DIR + "/bioinfo/Genomic/Sofware/GlimmerHMM/bin:" + self.config.SOFTWARE_DIR + "/bioinfo/annotation/fastree:" + self.config.SOFTWARE_DIR + "/bioinfo/align/ncbi-blast-2.3.0+/bin:" + self.config.SOFTWARE_DIR + "/bioinfo/align/clustalw-2.1/src:" + self.config.SOFTWARE_DIR + "/bioinfo/phylogenetic/muscle-3.8.31-release"
        self.similarity_file = os.path.join(self.config.SOFTWARE_DIR, "bioinfo/Genomic/Sofware/antismash-5.1.2/similarity.txt")
        self.lib = self.config.SOFTWARE_DIR + "/program/Python35/lib"
        self.set_environ(PATH=self.path, LD_LIBRARY_PATH=self.lib)
        self.set_environ(
            PATH=self.config.SOFTWARE_DIR + "/bioinfo/align/meme/meme/bin:" + self.config.SOFTWARE_DIR + "/mnt/ilustre/users/sanger-dev/meme/libexec/meme-4.11.2")

    def run_antismash(self):
        """
        运行run_antismash.py进行次级代谢分析注释
        :return:
        """
        if self.option("genome_seq").is_set:
            input_file = self.option("genome_seq").prop['path']
            self.file_name = os.path.splitext(os.path.basename(self.option("genome_seq").prop['path']))[0]
        elif self.option("genome_gbk").is_set:
            input_file = self.option("genome_gbk").prop['path']
            self.file_name = os.path.splitext(os.path.basename(self.option("genome_gbk").prop['path']))[0]
        if os.path.exists(self.work_dir + "/" + self.file_name):
            shutil.rmtree(self.work_dir + "/" + self.file_name)
        cmd = '{}antismash -c 10 --taxon {} {} --output-dir {}  --hmmdetection-strictness {} --cb-knownclusters --genefinding-tool {}'.format(
            "program/Python35/bin/", self.option("taxon"), input_file, self.file_name, self.option("strictness"),self.option("predict"))
        self.logger.info("运行antismash进行次级代谢分析注释")
        command = self.add_command("run_antismash", cmd, ignore_error=True).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_antismash运行完成")
            self.run_get_information()
        elif command.return_code == 1:
            self.logger.info("run_antismash运行完成,无结果")
        else:
            self.set_error("run_antismash运行出错!", code="31200101")

    def run_get_information(self):
        if os.path.exists(self.work_dir + "/" + self.file_name + "/knownclusterblast"):
            self.knowncluster()
            file_list = []
            for file in os.listdir(self.work_dir + "/" + self.file_name):
                if len(file.split(".")) > 2:
                    if file.split(".")[-2].startswith("region"):
                        file_list.append(file)
            file_list.sort()
            if file_list[0].startswith("Scaffold"):
                if file_list[0].split(".")[0] == "Scaffold":
                    pass
                else:
                    file_list.sort(key = lambda x:(int(x.split(".")[0].replace("Scaffold",""))))
            elif file_list[0].startswith("Chromosome"):
                if file_list[0].split(".")[0] == "Chromosome":
                    pass
                else:
                    file_list.sort(key=lambda x: (int(x.split(".")[0].replace("Chromosome", ""))))
            else:
                if file_list[0].split(".")[0] == "Plasmid":
                    pass
                else:
                    file_list.sort(key=lambda x: (str(x.split(".")[0].replace("Plasmid", ""))))
            if re.match("^[Cc]hromosome.*", self.file_name) or re.match("^[Pp]lasmid.*", self.file_name):
                file_name = self.file_name[0].lower() + self.file_name[1:] + "_antismash_anno.xls"
            else:
                file_name = "antismash_anno.xls"
            tidy_file = os.path.join(self.output_dir, file_name)
            gene_file = os.path.join(self.output_dir, "gene_antismash.xls")
            with open(tidy_file, "w") as f, open(gene_file, "w") as g:
                f.write("Cluster ID\tType\tStart\tEnd\tGene No.\tGenes\tpredicted_structure\tMost Similar Cluster\tSimilarity\tMIBiG accession\n")
                g.write("Gene ID\tType\tMost Similar Cluster\tLocation\tCluster ID\n")
                for i in range(len(file_list)):
                    a =  open(self.work_dir + "/" + self.file_name + "/" + file_list[i])
                    b = a.read()
                    a.close()
                    # cluster name
                    cluster = "cluster" + str(i + 1)
                    locus = file_list[i].split(".")[0]
                    locus1 = locus
                    #link gbk
                    if os.path.exists(self.output_dir + '/gbk/'):
                        pass
                    else:
                        os.mkdir(self.output_dir + '/gbk/')
                    if os.path.exists(self.output_dir + '/gbk/' +locus + "_" + cluster):
                        os.remove(self.output_dir + '/gbk/' +locus + "_" + cluster)
                    os.link(self.work_dir + "/" + self.file_name + "/" + file_list[i],self.output_dir + "/gbk/" + locus + "_" + cluster)
                    gene_num = 0
                    gene_name_list = []
                    c = b.split("\n")
                    # 提取Most Similar Cluster
                    if file_list[i].split(".")[1].startswith("region00"):
                        cluster_name = locus + "_c" + file_list[i].split(".")[1].replace("region00","")
                    elif file_list[i].split(".")[1].startswith("region0"):
                        cluster_name = locus + "_c" + file_list[i].split(".")[1].replace("region0","")
                    else:
                        cluster_name = locus + "_c" + file_list[i].split(".")[1].replace("region", "")
                    similarity = self.knowncluster_s[cluster_name]
                    msc = self.knowncluster_c[cluster_name]
                    accession = self.knowncluster_b[cluster_name]
                    # 提取type start end
                    cluster_type = b.partition(" region ")[2].partition("/product=\"")[2].partition("\"")[0]
                    if "\n" in cluster_type:
                        cluster_type = re.sub(r"\n +"," ",cluster_type)
                    start = b.partition("Orig. start  :: ")[2].partition("\n")[0]
                    end = b.partition("Orig. end    :: ")[2].partition("\n")[0]
                    #提取 产物结构
                    if b.partition("SMILES=")[1]:
                        core_structure = b.partition("SMILES=")[2].partition("\n")[0].replace("\"", "")
                    else:
                        core_structure = "-"
                    if core_structure == "" or core_structure == " ":
                        core_structure = "-"
                    # 提取 gene 信息
                    d = b.split(" gene ")
                    for j in d:
                        if j.partition("locus_tag=")[1]:
                            gene_name = j.partition("locus_tag=")[2].partition("\n")[0].replace("\"", "")
                            gene_num += 1
                            gene_name_list.append(gene_name)
                            g.write(gene_name + "\t" + cluster_type + "\t" + msc + "\t" + locus1 + "\t" + cluster  + "\n")
                    if str(start) == "0":
                        start = "1"
                    f.write(locus1 + "_" + cluster + "\t" + cluster_type + "\t" + str(start) + "\t" + str(end) + "\t" + str(gene_num) + "\t" + ";".join(gene_name_list) +
                            "\t" + core_structure + "\t" + msc + "\t" + str(similarity) + "\t" + accession + "\n")
                    gene_num = 0

    def knowncluster(self):
        similarity_data = {}
        with open(self.similarity_file) as f:
            a = f.readlines()
            for i in a:
                if i:
                    similarity_data[i.strip("\r\n").split("\t")[0]] = i.strip("\r\n").split("\t")[1]
        self.knowncluster_s = {}
        self.knowncluster_c = {}
        self.knowncluster_b = {}
        knowncluster_dir = os.path.join(self.work_dir, self.file_name + "/knownclusterblast")
        allclusterfiles = os.listdir(knowncluster_dir)
        for each in allclusterfiles:
            clusterfile = os.path.join(knowncluster_dir, each)
            self.logger.info(clusterfile)
            if os.path.isfile(clusterfile):
                cluster_name = each.replace(".txt", "") # contig115_c1
                # 原软件中提取获得相似性和knowncluster核心方法
                clusterblastfile = open(clusterfile, "r")
                clusterblastfile = clusterblastfile.read()
                clusterblastfile = clusterblastfile.replace("\r", "\n")
                hitlines = ((clusterblastfile.split("Significant hits: \n")[1]).split("\nDetails:")[0]).split("\n")
                firsthit = hitlines[0]  # 只获取第一条
                if firsthit != "":
                    accession = firsthit.partition(" ")[2].partition("\t")[0]
                    self.knowncluster_b[cluster_name] = accession
                    hitdetails = [details for details in (clusterblastfile.split("\nDetails:")[1]).split(">>")[1:] if accession in details][0]
                    simila_cluster = hitdetails.partition("Source: ")[2].partition("\n")[0]
                    simila_type = hitdetails.partition("Type: ")[2].partition("\n")[0]
                    # Get total number of genes in hitcluster
                    nr_hitclustergenes = len(
                        clusterblastfile.partition("of subject cluster:\n")[2].partition("\n\n")[0].split("\n"))
                    # Get percentage
                    blasthits = hitdetails.partition(" %coverage, e-value):\n")[2].partition("\n\n")[0].split("\n")
                    nr_genes_with_hits = len(set([hit.split("\t")[1] for hit in blasthits]))
                    percentage_genes_with_hits = int((float(nr_genes_with_hits) / float(similarity_data[accession])) * 100)
                    self.knowncluster_s[cluster_name] = percentage_genes_with_hits
                    self.knowncluster_c[cluster_name] = simila_cluster
                else:
                    self.knowncluster_s[cluster_name] = "-"
                    self.knowncluster_c[cluster_name] = "-"
                    self.knowncluster_b[cluster_name] = "-"

    def run(self):
        super(AntismashV5Tool, self).run()
        self.run_antismash()
        #self.run_get_information()
        self.set_output()

    def set_output(self):
        files = os.listdir(self.output_dir)
        for file in files:
            if re.search(r'antismash_anno.xls$', file):
                with open(self.output_dir + '/' + file) as p:
                    lines = p.readlines()
                    if len(lines) == 1:
                        os.remove(self.output_dir + '/' + file)
                    else:
                        pass
        if os.path.exists(self.output_dir + "/" + self.file_name + ".zip"):
            os.remove(self.output_dir + "/" + self.file_name + ".zip")
        if os.path.exists(self.output_dir + "/" + "gene_antismash.xls"):
            os.link(self.work_dir + "/" + self.file_name + "/" + self.file_name + ".zip", self.output_dir + "/" + self.file_name + ".zip")
        self.end()
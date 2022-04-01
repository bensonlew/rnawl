# -*- coding: utf-8 -*-
# __author__ = 'zzg'
# last_modify: 2020.08.27
from biocluster.agent import Agent
from biocluster.tool import Tool
import os, re,shutil
from biocluster.core.exceptions import OptionError
import subprocess
import unittest
import pandas as pd


class AntismashAgent(Agent):
    """
    功能：antismash注释,即次级代谢分析
    输入：组装得到的基因组序列(必须是*.fasta后缀)或者gbk文件
    比对软件：antismash-4.0.2
    """

    def __init__(self, parent):
        super(AntismashAgent, self).__init__(parent)
        options = [
            {"name": "input_file", "type": "infile", "format": "tool_lab.no_empty"},
            {"name": "input_format", "type": "string", "default": "gbk"},
            {"name": "taxon", "type": "string", "default": "bacteria"}, # 选择物种为真菌fungi时不能输入fasta文件
            {"name": "result_stat", "type": "outfile", "format": "sequence.profile_table"}
        ]
        self.add_option(options)
        self._memory_increase_step = 50  # 每次重运行增加内存50G

    def check_options(self):
        if self.option("input_format") in ["fasta"]:
            if  self.option("taxon") in ["fungi"]:
                raise OptionError("选择物种为真菌fungi时不能输入fasta文件", code="42200301")
        if self.option("input_format") in ["gbk"]:
            with open(self.option("input_file").prop['path'], 'rb') as r:
                file = r.readlines()
                locus = re.match("^LOCUS\s+(\S+)\s+(\w+) bp\s+\w+\s+(\w+).*$", file[0])
                locus_nu = 0
                gene_nu = 0
                genome_name = ""
                genome_length = ""
                genome_type = ""
                feature = ""
                if locus:
                    genome_name = locus.group(1)
                    genome_length = locus.group(2)
                    genome_type = locus.group(3)
                    locus_nu = 1
                else:
                    raise OptionError("文件格式错误，Genbank文件应该含有具体LOCUS信息", code="42200302")
                for f in file[1:]:
                    if re.match("^FEATURES\s+.*$", f):
                        feature = "ok"
                    if re.match("^ {5}CDS.*$", f):
                        gene_nu = gene_nu + 1
                    locus = re.match("^LOCUS\s+(\w+)\s+(\w+) bp\s+\w+\s+(\w+).*$", f)
                    if locus:
                        locus_nu = locus_nu + 1
                        genome_name = genome_name + "," + locus.group(1)
                        genome_length = str(genome_length) + "," + str(locus.group(2))
                        genome_type = genome_type + "," + locus.group(3)
                        if feature == "":
                            raise OptionError("文件格式错误，缺少FEATURES信息")
                        gene_number = ""
                        gene_number = gene_number + "\t" + str(gene_nu)
                        gene_nu = 0
        if self.option("input_format") in ["fasta"]:
            file_path = self.option("input_file").prop['path']
            if file_path.split(".")[-1] != "fasta" and file_path.split(".")[-1] != "Fasta" and file_path.split(".")[-1] != "FASTA":
                raise OptionError("文件格式错误", code="44000501")
        return True

    def set_resource(self):
        self._cpu = 10
        self._memory = '150G'    

    def end(self):
        super(AntismashAgent, self).end()


class AntismashTool(Tool):
    def __init__(self, config):
        super(AntismashTool, self).__init__(config)
        self.file_name = ""
        self.python = "/program/Python/bin/"
        self.antismash = os.path.join(self.config.SOFTWARE_DIR, "bioinfo/Genomic/Sofware/antismash-4.0.2/")
        self.set_environ(
            PATH=self.config.SOFTWARE_DIR + '/program/Python/bin/:' + self.config.SOFTWARE_DIR + '/bioinfo/align/hmmer-3.1b2-linux-intel-x86_64/binaries:' + self.config.SOFTWARE_DIR + '/bioinfo/Genomic/Sofware/glimmer3.02/bin:' + self.config.SOFTWARE_DIR + '/bioinfo/Genomic/Sofware/smrtanalysis_2.3.0/install/smrtanalysis_2.3.0.140936/analysis/bin/:' + self.config.SOFTWARE_DIR + '/bioinfo/Genomic/Sofware/Islander_software/bin:' + self.config.SOFTWARE_DIR + '/bioinfo/Genomic/Sofware/hmmer-2.3.2/bin:' + self.config.SOFTWARE_DIR + '/bioinfo/Genomic/Sofware/GlimmerHMM/bin:' + self.config.SOFTWARE_DIR + '/bioinfo/align/mafft-7.299-with-extensions/bin:' + self.config.SOFTWARE_DIR + '/bioinfo/align/ncbi-blast-2.3.0+/bin:' + self.config.SOFTWARE_DIR + '/bioinfo/align/clustalw-2.1/src:' + self.config.SOFTWARE_DIR + '/program/sun_jdk1.8.0/bin/:' + self.config.SOFTWARE_DIR +  "/bioinfo/phylogenetic/muscle-3.8.31-release")
        self.set_environ(LD_LIBRARY_PATH=self.config.SOFTWARE_DIR + '/program/Python/lib')

    def run_antismash(self):
        """
        运行run_antismash.py进行次级代谢分析注释
        :return:
        """
        input_file = self.option("input_file").prop['path']
        self.file_name = os.path.splitext(os.path.basename(self.option("input_file").prop['path']))[0]
        cmd = '{}python {}run_antismash.py -c 10 --taxon {} {}'.format(self.python, self.antismash,
                                                                       self.option("taxon"), input_file)
        #cmd += ' --knownclusterblast --transatpks_da --nclusters 1 '
        cmd += ' --knownclusterblast --nclusters 1 --smcogs '
        cmd += ' --disable-embl --disable-genbank ' #减少文件输出
        self.logger.info("运行run_antismash.py进行次级代谢分析注释")
        command = self.add_command("run_antismash", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_antismash运行完成")
        else:
            self.set_error("run_antismash运行出错!", code="31200101")

    def run_get_information(self):
        """
        提取次级代谢注释信息，用于页面展现
        """
        if os.listdir(self.work_dir + "/" + self.file_name + "/txt/"):
            self.knowncluster()
            antismash_files = self.work_dir + "/" + self.file_name + "/txt/" + "*_BGC.txt"
            antismash_path = self.work_dir + "/BGC.txt"
            os.system('cat ' + antismash_files + ' >' + antismash_path)
            if os.path.exists(antismash_path):
                '''
                if re.match("^[Cc]hromosome.*", self.file_name) or re.match("^[Pp]lasmid[A-Z].*", self.file_name):
                    file_name = self.file_name[0].lower() + self.file_name[1:] + "_antismash_anno.xls"
                else:
                '''
                file_name = "antismash_anno.xls"
                tidy_file = os.path.join(self.output_dir, file_name)
                with open(antismash_path, "r") as a, open(tidy_file, "w") as t:
                    sort_new = {}
                    id = []
                    t.write("Cluster ID\tType\tStart\tEnd\tGene No.\tGenes\tpredicted_structure\tMost Similar Cluster\tSimilarity\tsignature_genes\n")
                    anti = a.readlines()[1:]
                    for line in anti:
                        if line.startswith("BGC ID"):
                            pass
                        else:
                            arr = line.strip().split("\t")
                            if re.match("no_tag_f_[0-9]*", arr[4])  or re.match("cluster_[0-9]*_[a-zA-Z]*[0-9]*", arr[4]):  #or re.match("ctg[0-9]*_[0-9]*", arr[4])
                                pass
                            else:
                                cluster_id = arr[0].rsplit("_c")[0] + "_Cluster" + arr[0].rsplit("_c")[1]
                                cluster_id_tmp = re.subn("^[Cc]hromosome", "Chr", cluster_id)[0]
                                cluster_id = re.subn("^[Pp]lasmid", "p", cluster_id_tmp)[0]
                                cluster_type = arr[1].replace(";", "-")
                                start = int(arr[3].split(";")[0]) + 1
                                end = arr[3].split(";")[1]
                                arr4 = re.subn(";no_tag_f_[0-9]*", ";", arr[4])[0]
                                arr[4] = re.subn("^no_tag_f_[0-9]*;", "", arr4)[0]
                                '''
                                arr4 = re.subn(";ctg[0-9]*_[0-9]*", ";", arr[4])[0]
                                arr[4] = re.subn("^ctg[0-9]*_[0-9]*;", "", arr4)[0]
                                '''
                                arr4 = re.subn(";cluster_[0-9]*_[a-zA-Z]*[0-9]*", "", arr4)[0]
                                arr[4] = re.subn("^cluster_[0-9]*_[a-zA-Z]*[0-9]*;", "", arr4)[0]
                                gene_nu = len(arr[4].split(";"))
                                gene = arr[4]
                                signature_genes = arr[7]
                                # gene = arr4.replace("_1;", ";").rstrip("_1")
                                cluster_ori = "genecluster" + arr[0].rsplit("_c")[1]
                                cluster_name = 'cluster'+arr[0].rsplit("_c")[1]
                                predicted_structure = cluster_name + ".png"
                                self.logger.info(">>>>>>>>>>>>")
                                png_dir = self.work_dir + "/" + self.file_name + "/structures/"
                                self.logger.info(png_dir + predicted_structure)
                                if not os.path.exists(png_dir + cluster_ori+'.png'):
                                    predicted_structure = "-"
                                msc = self.knowncluster_c[cluster_name]
                                similarity = self.knowncluster_s[cluster_name]
                                new = cluster_id + "\t" + cluster_type + "\t" + str(start) + "\t" + end + "\t" + str(
                                    gene_nu) + "\t" + gene + "\t" + predicted_structure + "\t"+ msc + "\t" + str(similarity) + "\t" + signature_genes + "\n"
                                idn = int(arr[0].rsplit("_c")[1])
                                sort_new[idn] = new
                                id.append(idn)
                    for i in sorted(id):
                        t.write(sort_new[i])

                with open(tidy_file) as t:
                    lines = t.readlines()
                if len(lines) < 2:
                #if True:
                    #self.option('result_stat', tidy_file)
                    os.remove(tidy_file)
                    self.logger.info('*_antismash_anno.xls 结果为空')

                os.system(
                    "cd " + self.work_dir + ";tar czvf antismash.tar.gz " + self.file_name + ";mv antismash.tar.gz output/")
                self.option('result_stat', tidy_file)
            self.creat_gene_table(tidy_file)

    def knowncluster(self):
        self.knowncluster_s = {}
        self.knowncluster_c = {}
        knowncluster_dir = os.path.join(self.work_dir, self.file_name +"/knownclusterblast")
        allclusterfiles = os.listdir(knowncluster_dir)
        for each in allclusterfiles:
            clusterfile = os.path.join(knowncluster_dir, each)
            self.logger.info(clusterfile)
            if os.path.isfile(clusterfile):
                cluster_name = each.replace(".txt","")
                # 原软件中提取获得相似性和knowncluster核心方法
                clusterblastfile = open(clusterfile,"r")
                clusterblastfile = clusterblastfile.read()
                clusterblastfile = clusterblastfile.replace("\r","\n")
                hitlines = ((clusterblastfile.split("Significant hits: \n")[1]).split("\nDetails:")[0]).split("\n")
                firsthit= hitlines[0] # 只获取第一条
                if firsthit !="":
                    accession = firsthit.partition(" ")[2].partition("\t")[0]
                    hitdetails = [details for details in (clusterblastfile.split("\nDetails:")[1]).split(">>")[1:] if accession in details][0]
                    simila_cluster = hitdetails.partition("Source: ")[2].partition("\n")[0]
                    simila_type = hitdetails.partition("Type: ")[2].partition("\n")[0]
                    #Get total number of genes in hitcluster
                    nr_hitclustergenes = len(clusterblastfile.partition("of subject cluster:\n")[2].partition("\n\n")[0].split("\n"))
                    #Get percentage
                    blasthits = hitdetails.partition(" %coverage, e-value):\n")[2].partition("\n\n")[0].split("\n")
                    nr_genes_with_hits = len(set([hit.split("\t")[1] for hit in blasthits]))
                    percentage_genes_with_hits = int((float(nr_genes_with_hits) / float(nr_hitclustergenes)) * 100)
                    self.knowncluster_s[cluster_name] = percentage_genes_with_hits
                    self.knowncluster_c[cluster_name] = simila_cluster + "," + simila_type
                else:
                    self.knowncluster_s[cluster_name] = "-"
                    self.knowncluster_c[cluster_name] = "-"

    def creat_gene_table(self, tidy_file):
        tidy_table = pd.read_table(tidy_file, header=0, sep="\t")
        stat_table = tidy_table[["Genes","Cluster ID","Type","Most Similar Cluster"]]
        stat_table = stat_table.drop("Genes", axis=1).join(stat_table["Genes"].str.split(";", expand=True).stack().reset_index(\
            level=1, drop=True).rename("Gene ID"))
        split_table = stat_table["Cluster ID"].str.split("_Cluster", 1, expand=True)
        split_table[1] = "Cluster" + split_table[1]
        split_table.columns = ["Location","Cluster ID"]
        final_table = pd.concat([stat_table[["Gene ID","Type","Most Similar Cluster"]],split_table],axis=1)
        final_table.to_csv(self.output_dir + "/gene_antismash.xls",sep="\t",index=False)

    def get_gene_information(self):
        if os.listdir(self.work_dir + "/" + self.file_name + "/txt/"):
            gene_files = self.work_dir + "/" + self.file_name + "/txt/" + "*_gene.txt"
            gene_path = self.work_dir + "/gene.txt"
            os.system('cat ' + gene_files + ' >' + gene_path)
            if os.path.exists(gene_path):
                file_name = "gene_list.xls"
                gene_list_file = os.path.join(self.output_dir, file_name)
                with open(gene_path, "r") as y, open(gene_list_file, "w") as g:
                    f = open(self.work_dir + "/" + self.file_name + "/smcogs/" + "smcogs.txt")
                    a = f.read()
                    f.close()
                    b = a.split(">>")
                    g.write("Gene ID\tGene start\tGene end\tGene strand\tsmcog\ttype\n")
                    gene_info = y.readlines()[1:]
                    for line in gene_info:
                        if line.startswith("gene ID"):
                            pass
                        else:
                            arr = line.strip().split("\t")
                            gene_name = arr[0]
                            start = str(arr[1])
                            end = str(arr[2])
                            strand = str(arr[3])
                            if re.search(arr[0]+"\n", a):
                                for i in range(len(b)):
                                    if re.search(arr[0]+"\n", b[i]):
                                        if re.search("transport", b[i]):
                                            c = b[i].split()
                                            g.write(gene_name + "\t" + start + "\t" + end + "\t" + strand + "\t" + c[10].split(",")[0] + "\t" + "transport" + "\n")
                                        elif re.search("regulator", b[i]):
                                            d = b[i].split()
                                            g.write(gene_name + "\t" + start + "\t" + end + "\t" + strand + "\t" + d[10].split(",")[0] + "\t" + "regulator" + "\n")
                                        elif (not re.search("transport", b[i])) and (not re.search("regulator", b[i])) and (re.search("smCOG", b[i])):
                                            e = b[i].split()
                                            g.write(gene_name + "\t" + start + "\t" + end + "\t" + strand + "\t" + e[10].split(",")[0] + "\t" + "additional" + "\n")
                                    else:
                                        continue
                            else:
                                g.write(gene_name + "\t" + start + "\t" + end + "\t" + strand + "\t" + "-" + "\t" + "-" + "\n")

    def run(self):
        super(AntismashTool, self).run()
        self.run_antismash()
        self.run_get_information()
        self.get_gene_information()
        self.set_output()
        self.end()

    def set_output(self):
        files =os.listdir(self.output_dir)
        for file in files:
            if re.search(r'antismash_anno.xls$',file):
                with open (self.output_dir + '/' + file) as p:
                    lines =p.readlines()
                    if len(lines) == 1:
                        os.remove(self.output_dir + '/' + file)
                    else:
                        pass
        dd =os.listdir(self.output_dir)
        if len(dd) == 1:
            for d in dd:
                os.remove(self.output_dir + '/' +d)
        else:
            pass
        '''
        新增核心产物结果图片结果
        '''
        self.logger.info("link core structure>>>>>>")
        png_dir = self.work_dir + "/" + self.file_name + "/structures/"
        if os.path.exists(png_dir):
            png_files = os.listdir(png_dir)
            new_dir = os.path.join(self.output_dir, "core_structures")
            if not os.path.exists(new_dir):
                os.mkdir(new_dir)
            for eachfile in png_files:
                if re.match("genecluster\d+.png", eachfile):
                    newname = eachfile.replace("genecluster","cluster")
                    old_file = os.path.join(png_dir, eachfile)
                    new_file = os.path.join(new_dir, newname)
                    if os.path.exists(new_file):
                        os.remove(new_file)
                    os.link(old_file, new_file)
# -*- coding: utf-8 -*-
# __author__ = 'zhujuan'
# last_modify: 2018.02.07
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
            {"name": "genome_seq", "type": "infile", "format": "sequence.fasta"},
            # 组装得到的基因组文件，尾缀必须是.fasta， 该步先进行基因预测（默认glimmer3）
            {"name": "genome_gbk", "type": "infile", "format": "gene_structure.gbk"},  # 样品基因组的gbk格式文件，无需进行预测
            {"name": "taxon", "type": "string", "default": "bacteria"},
            # 基因组序列来源：bacteria,fungi #目前软件中当输入为*.fasta，且taxon 为 fungi时，没有对应的基因预测方法
            {"name": "result_stat", "type": "outfile", "format": "sequence.profile_table"},
        ]
        self.add_option(options)
        #self.queue = 'SANGER'  # 去掉投递 by zzg 20201204
        self._memory_increase_step = 50  # 每次重运行增加内存50G by qingchen.zhang @ 20200916

    def check_options(self):
        if not self.option("genome_seq").is_set and not self.option("genome_gbk").is_set:
            raise OptionError("必须设置参数genome_seq或者genome_gbk", code="31200101")
        return True

    def set_resource(self):
        self._cpu = 10
        self._memory = '75G'

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
        if self.option("genome_seq").is_set:
            input_file = self.option("genome_seq").prop['path']
            self.file_name = os.path.splitext(os.path.basename(self.option("genome_seq").prop['path']))[0]
        elif self.option("genome_gbk").is_set:
            input_file = self.option("genome_gbk").prop['path']
            self.file_name = os.path.splitext(os.path.basename(self.option("genome_gbk").prop['path']))[0]
        cmd = '{}python {}run_antismash.py -c 16 --taxon {} {}'.format(self.python, self.antismash,
                                                                       self.option("taxon"), input_file)
        #cmd += ' --knownclusterblast --transatpks_da --nclusters 1 '
        cmd += ' --knownclusterblast --nclusters 1 '
        cmd += ' --disable-embl --disable-genbank ' #减少文件输出
        self.logger.info("运行run_antismash.py进行次级代谢分析注释")
        command = self.add_command("run_antismash", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_antismash运行完成")
            self.run_get_information()
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
                if re.match("^[Cc]hromosome.*", self.file_name) or re.match("^[Pp]lasmid[A-Z].*", self.file_name):
                    file_name = self.file_name[0].lower() + self.file_name[1:] + "_antismash_anno.xls"
                else:
                    file_name = "antismash_anno.xls"
                tidy_file = os.path.join(self.output_dir, file_name)
                with open(antismash_path, "r") as a, open(tidy_file, "w") as t:
                    sort_new = {}
                    id = []
                    t.write("Cluster ID\tType\tStart\tEnd\tGene No.\tGenes\tpredicted_structure\tMost Similar Cluster\tSimilarity\n")
                    anti = a.readlines()[1:]
                    for line in anti:
                        if line.startswith("BGC ID"):
                            pass
                        else:
                            arr = line.strip().split("\t")
                            if re.match("no_tag_f_[0-9]*", arr[4]) or re.match("ctg[0-9]*_[0-9]*", arr[4]) or re.match("cluster_[0-9]*_[a-zA-Z]*[0-9]*", arr[4]):
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
                                arr4 = re.subn(";ctg[0-9]*_[0-9]*", ";", arr[4])[0]
                                arr[4] = re.subn("^ctg[0-9]*_[0-9]*;", "", arr4)[0]
                                arr4 = re.subn(";cluster_[0-9]*_[a-zA-Z]*[0-9]*", "", arr4)[0]
                                arr[4] = re.subn("^cluster_[0-9]*_[a-zA-Z]*[0-9]*;", "", arr4)[0]
                                gene_nu = len(arr[4].split(";"))
                                gene = arr[4]
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
                                    gene_nu) + "\t" + gene + "\t" + predicted_structure + "\t"+ msc + "\t" + str(similarity) +"\n"
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
                    self.end()

                os.system(
                    "cd " + self.work_dir + ";tar czvf antismash.tar.gz " + self.file_name + ";mv antismash.tar.gz output/")
                self.option('result_stat', tidy_file)
            self.creat_gene_table(tidy_file)
        #        self.end()
        else:
            self.end()

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
        split_table = stat_table["Cluster ID"].str.split("_",1,expand=True)
        split_table.columns = ["Location","Cluster ID"]
        final_table = pd.concat([stat_table[["Gene ID","Type","Most Similar Cluster"]],split_table],axis=1)
        final_table.to_csv(self.output_dir + "/gene_antismash.xls",sep="\t",index=False)

    def run(self):
        super(AntismashTool, self).run()
        self.run_antismash()
        self.set_output()

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
        self.end()


class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            #"id": "CreatTable" + str(random.randint(1, 10000)),
            "id": "Antismash_1",
            "type": "tool",
            "name": "annotation.antismash",
            "instant": False,
            "options": dict(
                genome_gbk="/mnt/ilustre/users/sanger-dev/workspace/20190308/Bacgenome_tsg_33557/Bacgenome/BacGbk/output/analysis_gbk/all.antismash.gbk",
                #genome_gbk="/mnt/ilustre/users/sanger-dev/sg-users/yuanshaohua/measurement/antismash/Y16952.3.region001.gbk",
                taxon="bacteria",

            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
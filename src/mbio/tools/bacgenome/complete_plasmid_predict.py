# -*- coding: utf-8 -*-
# __author__ = 'hao.gao'
# last_modify: 2021.04.03

import os, time
import re, shutil
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from mbio.packages.bac_comp_genome.common_function import link_dir
from Bio import SeqIO
import pandas as pd
from pandas.core.frame import DataFrame

class CompletePlasmidPredictAgent(Agent):
    """
    完成图质粒预测小工具
    采用conda进行
    """
    def __init__(self, parent):
        super(CompletePlasmidPredictAgent, self).__init__(parent)
        options = [
            {"name": "fasta", "type": "infile", "format": "sequence.fasta"},# 序列文件
            {"name": "status", "type": "infile", "format": "sequence.profile_table"},  # 记录完成图的成环状态
            {"name": "sample_name", "type": "string"},# 样本名
            {"name": "threshold", "type": "float", "default": 0.7},
            {"name": "outfmt", "type": "int", "default": 6},
            {"name": "evalue", "type": "float", "default": 1e-5},
            {"name": "max_target_seqs", "type": "int", "default": 1},
            {"name": "scf", "type": "outfile", "format": "sequence.fasta"},
            {"name": "table", "type": "outfile", "format": "sequence.profile_table"}
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("fasta").is_set:
            raise OptionError("必须输入序列文件")
        if not self.option("status").is_set:
            raise OptionError("必须输入status文件！")

    def set_resource(self):
        self._cpu = 4
        self._memory = '20G'

    def end(self):
        super(CompletePlasmidPredictAgent, self).end()

class CompletePlasmidPredictTool(Tool):
    def __init__(self, config):
        super(CompletePlasmidPredictTool, self).__init__(config)
        self.conda = self.config.SOFTWARE_DIR + "/program/miniconda3/"
        self.set_environ(PATH=self.conda +"bin")
        self.path = self.config.SOFTWARE_DIR + "/program/Python35/bin:"
        self.lib = self.config.SOFTWARE_DIR + "/program/Python35/lib"
        self.set_environ(PATH=self.path, LD_LIBRARY_PATH=self.lib)
        self.diamond = "/bioinfo/align/diamond-0.9.11/diamond"
        self.shell_path = "/program/sh"
        self.database = self.config.SOFTWARE_DIR + "/database"
        self.python = "/program/Python35/bin/python"
        self.script = self.config.PACKAGE_DIR + '/bacgenome/'


    def creat_shell(self, input_file, threshold):
        """
        重写一个shell，将运行的命令直接写入文件里面去
        :return:
        """
        self.fix_sh = os.path.join(self.work_dir, "plas_flow.sh")
        with open(self.fix_sh, "w") as w:
            w.write("#! /bin/bash\n")
            w.write("source {}\n".format(self.config.SOFTWARE_DIR + "/program/miniconda3/etc/profile.d/conda.sh"))
            w.write("conda activate plasflow-1.1\n")
            w.write("PlasFlow.py --input {} --output plasflow_predictions.tsv --threshold {} > plasflow.log\n".format(input_file,threshold))
            w.write("conda deactivate\n")

    def run_predict(self):
        """
        质粒鉴定
        :return:
        """
        self.creat_shell(self.option("fasta").prop['path'], self.option("threshold"))
        cmd = "{} {}".format(self.shell_path, self.fix_sh)
        self.logger.info(cmd)
        command = self.add_command("run_predict", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_predict运行完成！")
        else:
            self.set_error("run_predict运行运行出错!")
        self.plas = os.path.getsize(self.work_dir + "/plasflow_predictions.tsv_plasmids.fasta")
        return self.plas

    def run_dnaA(self):
        """
        dnaA注释
        :return:
        """
        cmd = "{} blastx -d {} -q {} -o {}  --max-target-seqs 1".format(self.diamond, self.database + "/dnaA/dnaA",self.option("fasta").prop['path'], self.work_dir + "/dnaA.blast.xls")
        self.logger.info(cmd)
        command = self.add_command("run_dnaa", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_dnaA运行完成！")
        else:
            self.set_error("run_dnaA运行运行出错!")

    def run_stat(self):
        dna_dict = {}
        if os.path.getsize(self.work_dir + "/dnaA.blast.xls") >0:
            with open(self.work_dir + "/dnaA.blast.xls", 'r') as t:
                liness = t.readlines()
                for line in liness:
                    lin = line.strip().split("\t")
                    if float(lin[2]) >= 90:
                        dna_dict[lin[0]] = 1
                    else:
                        dna_dict[lin[0]] = 0
        with open(self.work_dir + "/plasflow_predictions.tsv", "r") as f, open(self.option('status').prop['path'], "r") as g,open(self.output_dir+"/" +self.option("sample_name")+ ".stat.xls", "w") as s:
            s.write("sample\tknown_plasmid\tchromosome\tnew_plasmid\n")
            status_dict = {}
            lines =g.readlines()
            chr_num =0
            big_num = 0
            pla_num = 0
            for line in lines:
                lin = line.strip().split("\t")
                status_dict[lin[0]] = lin[1]
            lines2 = f.readlines()
            plasmid = 0
            chromosome = 0
            unclassified = 0
            plasmid_list = []
            fasta_list = []
            for line in lines2[1:]:
                lin = line.strip().split("\t")
                type = lin[5].split(".")[0]
                num = 6 + int(lin[4])
                de = ""
                if int(lin[3]) >= 2000 :
                    if type == "plasmid":
                        if int(lin[3]) >= 1000000 and lin[2] not in dna_dict:
                            de = 'known-Plasmid'
                            plasmid += 1
                            fasta_list.append([lin[2], status_dict[lin[2]], lin[3], "plasmid", de])
                        elif int(lin[3]) >= 1000000 and lin[2] in dna_dict:
                            if dna_dict[lin[2]] == 1:
                                de = 'Chromosome'
                                chromosome += 1
                                fasta_list.append([lin[2], status_dict[lin[2]], lin[3], "chromosome", de])
                            elif dna_dict[lin[2]] == 0:
                                de = 'known-Plasmid'
                                plasmid += 1
                                fasta_list.append([lin[2], status_dict[lin[2]], lin[3], "plasmid", de])

                        if int(lin[3]) < 1000000:
                            de = 'known-Plasmid'
                            plasmid += 1
                            fasta_list.append([lin[2], status_dict[lin[2]], lin[3], "plasmid", de])
                    elif type == "chromosome":
                        if int(lin[3]) >= 1000000:
                            de = 'Chromosome'
                            chromosome += 1
                            fasta_list.append([lin[2], status_dict[lin[2]], lin[3], "chromosome", de])
                        else:
                            if lin[2] in dna_dict:
                                if dna_dict[lin[2]] == 1:
                                    de = 'Chromosome'
                                    chromosome += 1
                                    fasta_list.append([lin[2], status_dict[lin[2]], lin[3], "chromosome", de])
                                elif dna_dict[lin[2]] == 0:
                                    de = 'New-Plasmid'
                                    unclassified += 1
                                    fasta_list.append([lin[2], status_dict[lin[2]], lin[3], "plasmid", de])
                            else:
                                de = 'New-Plasmid'
                                unclassified += 1
                                fasta_list.append([lin[2], status_dict[lin[2]], lin[3], "plasmid", de])
                    elif type == "unclassified":
                        if int(lin[3]) >= 1000000:
                            if lin[2] in dna_dict:
                                if dna_dict[lin[2]] == 1:
                                    de = 'Chromosome'
                                    chromosome += 1
                                    fasta_list.append([lin[2], status_dict[lin[2]], lin[3], "chromosome", de])
                                elif dna_dict[lin[2]] == 0:
                                    de = 'New-Plasmid'
                                    unclassified += 1
                                    fasta_list.append([lin[2], status_dict[lin[2]], lin[3], "plasmid", de])
                            else:
                                de = 'New-Plasmid'
                                unclassified += 1
                                fasta_list.append([lin[2], status_dict[lin[2]], lin[3], "plasmid", de])
                        else:
                            de = 'New-Plasmid'
                            unclassified += 1
                            fasta_list.append([lin[2], status_dict[lin[2]], lin[3], "plasmid", de])
                    if de == "Chromosome" and status_dict[lin[2]] in ["circle","Circle"]:
                        big_num += 1
                    if de == "New-Plasmid" or de == 'known-Plasmid':
                        if status_dict[lin[2]] in ["linear", "Linear"] and int(lin[3]) >= 10000:
                            pla_num +=1
                else:
                    if type == "plasmid":
                        de = 'known-Plasmid'
                        plasmid += 1
                        fasta_list.append([lin[2], status_dict[lin[2]], lin[3], "plasmid", de])
                if de =="New-Plasmid" or de == 'known-Plasmid':
                    for iterator in SeqIO.parse(self.option("fasta").prop['path'], "fasta"):
                        if iterator.id == lin[2] and status_dict[lin[2]] in ["circle","Circle"]:
                            plasmid_list.append(iterator)
                if de =="Chromosome":
                    chr_num +=1
            SeqIO.write(plasmid_list, self.work_dir+"/plasmid.fasta", "fasta")
            data = DataFrame(fasta_list)
            data.columns = ['seq_id', "status", "length", "type", 'anno_type']
            data["length"] = data["length"].astype(int)
            data.sort_values(["type", "length"], inplace=True,ascending=[True, False])
            list2 =[]
            for i in ['plasmid', 'chromosome']:
                data1 = data[data["type"] == i]
                if i == "chromosome":
                    if data1.shape[0] == 1:
                        data1['location'] = ["Chromosome"]
                    else:
                        list1 = []
                        for j in range(1, data1.shape[0] + 1):
                            list1.append("Chromosome" + str(j))
                        data1['location'] = list1
                    list2.append(data1)
                elif i == "plasmid":
                    list1 = []
                    dd ={1:"A", 2:"B",3:"C",4:"D",5:"E",6:"F",7:"G",8:"H",9:"I",10:"J",11:"K",12:"L",13:"M",14:"N",15:"O",16:"P",17:"Q",18:"R",19:"S"}
                    for j in range(1, data1.shape[0] + 1):
                        list1.append("Plasmid" + str(dd[j]))
                    data1['location'] = list1
                    list2.append(data1)
            all_data = pd.concat(list2)
            all_data.to_csv(self.output_dir+"/all.stat.xls", sep='\t', header=True, index=False,
                            columns=['seq_id', "status", "length", "type", 'anno_type', "location"])
            ass_list = []
            all_data['length'] = all_data['length'].astype(int)
            all_data.sort_values(['length'], inplace=True, ascending=False)
            circle_list =[]
            for index, row in all_data.iterrows():
                for iterator in SeqIO.parse(self.option("fasta").prop['path'], "fasta"):
                    if iterator.id == row["seq_id"] and row["status"] == "circle":
                        circle_list.append(iterator.id)
                        iterator.id = row["location"]
                        ass_list.append(iterator)
            SeqIO.write(ass_list, self.output_dir+"/"+self.option('sample_name')+".assemble.fasta", "fasta")
            fail_list = []
            for iterator in SeqIO.parse(self.option("fasta").prop['path'], "fasta"):
                if iterator.id not in circle_list:
                    fail_list.append(iterator)
            if len(fail_list) >0:
                SeqIO.write(fail_list, self.output_dir + "/.failing.fasta", "fasta")
            s.write("{}\t{}\t{}\t{}\n".format(self.option('sample_name'), plasmid, chromosome, unclassified))
            return chr_num,big_num,pla_num


    def run_anno(self):
        """
        质粒注释
        :return:
        """
        cmd = "{} -query {} -db {} -out plasmid_blast.txt -outfmt {} -evalue {} -max_target_seqs {} -max_hsps 1".format(
            "bioinfo/align/ncbi-blast-2.3.0+/bin/blastn", self.work_dir+"/plasmid.fasta", self.database + "/PLSDB/plsdb.fna",
            self.option("outfmt"), self.option("evalue"), self.option("max_target_seqs"))
        self.logger.info(cmd)
        command = self.add_command("run_anno", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_anno运行完成！")
        else:
            self.set_error("run_anno运行运行出错!")

    def run_anno_stat(self):
        cmd = "{} {}plasmid_anno.py -i {} -d {} -f {} -t {} -o {}".format(self.python, self.script, self.work_dir+"/plasmid_blast.txt", self.database + "/PLSDB/plsdb.tsv", self.option("fasta").prop['path'], self.output_dir+"/all.stat.xls", self.output_dir + "/"+self.option("sample_name") +".plasmid_anno.xls")
        self.logger.info(cmd)
        command = self.add_command("run_anno_stat", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_anno_stat运行完成！")
        else:
            self.set_error("run_anno_stat运行运行出错!")

    def set_output(self):
        """
        设置结果文件目录
        """
        self.option("scf", self.output_dir+"/" + self.option('sample_name') +".assemble.fasta")
        self.option("table", self.output_dir+"/all.stat.xls")
        self.logger.info("生成结果文件完成")

    def get_num(self):
        list =[]
        for iterator in SeqIO.parse(self.option("fasta").prop['path'], "fasta"):
            list.append(iterator.id)
        return len(list)

    def run(self):
        """
        运行
        """
        super(CompletePlasmidPredictTool, self).run()
        self.run_predict()
        num = self.get_num()
        if num < 20:
            self.run_dnaA()
            chr_num, big_num, pla_num = self.run_stat()
            if chr_num == 0:
                self.set_error("样品{},染色体数 = 0 !".format(self.option("sample_name")))
            if big_num == 0:
                self.set_error("样品{},长序列(>2K)成环染色体数  = 0 !".format(self.option("sample_name")))
            if pla_num > 0:
                self.set_error("样品{},长序列(>2K)线性质粒 >0 !".format(self.option("sample_name")))
            if os.path.getsize(self.work_dir+"/plasmid.fasta") >0:
                self.run_anno()
                self.run_anno_stat()
            else:
                self.end()
            self.set_output()
            self.end()
        else:
            self.set_error("样品{},组装序列总数>20条，请核查原始数据质量!".format(self.option("sample_name")))


# -*- coding: utf-8 -*-
# __author__ = 'hao.gao'
# last_modify: 2021.04.03

import os, time
import re, shutil
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from mbio.packages.bac_comp_genome.common_function import link_dir
from pandas.core.frame import DataFrame
from Bio import SeqIO

class DraftPlasmidPredictAgent(Agent):
    """
    扫描图质粒预测小工具
    采用conda进行
    """
    def __init__(self, parent):
        super(DraftPlasmidPredictAgent, self).__init__(parent)
        options = [
            {"name": "fasta", "type": "infile", "format": "sequence.fasta"},# 序列文件
            {"name": "sample_name", "type": "string"},# 样本名
            {"name": "threshold", "type": "float", "default": 0.7},
            {"name": "outfmt", "type": "int", "default": 6},
            {"name": "evalue", "type": "float", "default": 1e-5},
            {"name": "max_target_seqs", "type": "int", "default": 1},
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("fasta").is_set:
            raise OptionError("必须输入序列文件")
        if not self.option("fasta").is_set:
            raise OptionError("必须输入样品名称")

    def set_resource(self):
        self._cpu = 4
        self._memory = '20G'

    def end(self):
        super(DraftPlasmidPredictAgent, self).end()

class DraftPlasmidPredictTool(Tool):
    def __init__(self, config):
        super(DraftPlasmidPredictTool, self).__init__(config)
        self.conda = self.config.SOFTWARE_DIR + "/program/miniconda3/"
        self.set_environ(PATH=self.conda +"bin")
        self.path = self.config.SOFTWARE_DIR + "/program/Python35/bin:"
        self.lib = self.config.SOFTWARE_DIR + "/program/Python35/lib"
        self.set_environ(PATH=self.path, LD_LIBRARY_PATH=self.lib)
        self.shell_path = "/program/sh"
        self.database = self.config.SOFTWARE_DIR + "/database"
        self.python = "/miniconda2/bin/python"
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

    def run_anno(self):
        """
        质粒注释
        :return:
        """
        cmd = "{} -query {}/plasflow_predictions.tsv_plasmids.fasta -db {} -out plasmid_blast.txt -outfmt {} -evalue {} -max_target_seqs {}".format(
            "bioinfo/align/ncbi-blast-2.3.0+/bin/blastn", self.work_dir, self.database + "/PLSDB/plsdb.fna",
            self.option("outfmt"), self.option("evalue"), self.option("max_target_seqs"))
        self.logger.info(cmd)
        command = self.add_command("run_anno", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_anno运行完成！")
        else:
            self.set_error("run_anno运行运行出错!")

    def run_anno_stat(self):
        cmd = "{} {}plasmid_anno.py -i {} -d {} -f {} -t {} -o {}".format(self.python, self.script, self.work_dir+"/plasmid_blast.txt", self.database + "/PLSDB/plsdb.tsv", self.option("fasta").prop['path'], self.output_dir + "/all.stat.xls", self.output_dir + "/"+self.option("sample_name") +".plasmid_anno.xls")
        self.logger.info(cmd)
        command = self.add_command("run_anno_stat", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_anno_stat运行完成！")
        else:
            self.set_error("run_anno_stat运行运行出错!")

    def run_stat(self):
        with open(self.work_dir + "/plasflow_predictions.tsv", "r") as f, open(self.output_dir+"/" +self.option("sample_name")+ ".stat.xls", "w") as g, open(self.output_dir+"/" +self.option("sample_name")+ ".detail.xls", "w") as k:
            k.write("sample\tlocation\ttype\tphylum\tprobability\n")
            g.write("sample\tknown_plasmid\tchromosome\tnew_plasmid\n")
            lines = f.readlines()
            plasmid = 0
            chromosome = 0
            unclassified = 0
            list1 = []
            for line in lines[1:]:
                lin = line.strip().split("\t")
                num = 6+int(lin[4])
                if lin[5].split(".")[0] == "plasmid":
                    de = 'known-Plasmid'
                    plasmid += 1
                elif lin[5].split(".")[0] == "chromosome":
                    de = 'Chromosome'
                    chromosome += 1
                elif lin[5].split(".")[0] == "unclassified":
                    de = 'New-Plasmid'
                    unclassified += 1
                list1.append([lin[2], "-",lin[3], "-", de])
                k.write("{}\t{}\t{}\t{}\t{}\n".format(self.option('sample_name'), lin[2], de, lin[5].split(".")[1], lin[num]))
            g.write("{}\t{}\t{}\t{}\n".format(self.option('sample_name'), plasmid, chromosome, unclassified))
            data = DataFrame(list1)
            data.columns = ['seq_id', "satus", "length", "type", 'anno_type']
            data["length"] = data["length"].astype(int)
            data.sort_values(["length"], inplace=True, ascending= False)
            num = 0
            ass_list =[]
            list2 = []
            for index, row in data.iterrows():
                num +=1
                for iterator in SeqIO.parse(self.option("fasta").prop['path'], "fasta"):
                    if iterator.id == row["seq_id"]:
                        iterator.id = "Scaffold"+str(num)
                        ass_list.append(iterator)
                        list2.append("Scaffold"+str(num))
            SeqIO.write(ass_list, self.output_dir+"/"+self.option('sample_name')+".assemble.fasta", "fasta")
            data['location'] =list2
            data.to_csv(self.output_dir + "/all.stat.xls", sep='\t', header=True, index=False,
                            columns=['seq_id', "satus", "length", "type", 'anno_type', "location"])


    def set_output(self):
        """
        设置结果文件目录
        """
        self.logger.info("生成结果文件完成")

    def run(self):
        """
        运行
        """
        super(DraftPlasmidPredictTool, self).run()
        plas = self.run_predict()
        if plas == 0:
            self.run_stat()
            self.end()
        else:
            self.run_stat()
            self.run_anno()
            self.run_anno_stat()
        self.set_output()
        self.end()
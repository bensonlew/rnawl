# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
# version 1.0
# last_modify: 2021.04.19

import os, time
import re, shutil
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from mbio.packages.bac_comp_genome.common_function import link_dir
import json

class PhigaroAgent(Agent):
    """
    单个基因组phigaro预测前噬菌体
    """
    def __init__(self, parent):
        super(PhigaroAgent, self).__init__(parent)
        options = [
            {"name": "fasta", "type": "infile", "format": "sequence.fasta"},# 序列文件
            {"name": "sample_name", "type": "string"},# 样本名
            #{"name": "anno", "type": "infile", "format": "sequence.profile_table"},  # 注释总览表
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("fasta").is_set:
            raise OptionError("必须输入序列文件")

    def set_resource(self):
        self._cpu = 8
        self._memory = '20G'

    def end(self):
        super(PhigaroAgent, self).end()

class PhigaroTool(Tool):
    def __init__(self, config):
        super(PhigaroTool, self).__init__(config)
        self.path = self.config.SOFTWARE_DIR + "/program/Python35/bin"
        self.lib = self.config.SOFTWARE_DIR + "/program/Python35/lib"
        self.set_environ(PATH=self.path, LD_LIBRARY_PATH=self.lib)
        self.db = self.config.SOFTWARE_DIR + "/bioinfo/Genomic/Sofware/phigaro/config.yml"
        self.python = "/program/Python/bin/python"
        self.python2 = "/program/Python35/bin/phigaro"
        self.script = self.config.PACKAGE_DIR + '/bacgenome/'

    def run_phigaro(self):
        """
        运行phigaro
        :return:
        """
        cmd = "{} -f {} -p -e html tsv gff -o {} -c {} -d -t 4".format(self.python2, self.option("fasta").prop['path'], self.work_dir+"/phigaro", self.db)
        self.logger.info(cmd)
        command = self.add_command("run_phigaro", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_phigaro运行完成！")
        else:
            self.set_error("run_phigaro运行出错!")
        if os.path.exists(self.work_dir+"/phigaro"):
            self.phigaro_dir = self.work_dir+"/phigaro"
        else:
            self.phigaro_dir = self.work_dir
        for file in os.listdir(self.phigaro_dir):
            if file.endswith(".tsv"):
                self.phigaro_tsv = self.phigaro_dir + "/" + file
            if file.endswith(".gff3"):
                self.phigaro_gff3 = self.phigaro_dir + "/" + file
        with open(self.phigaro_tsv, "r") as f:
            lines = f.readlines()
            return len(lines)

    def run_stat(self):
        """
        运行run_stat,主要统计信息
        :return:
        """
        if os.path.exists(self.output_dir + "/" +self.option("sample_name")):
            shutil.rmtree(self.output_dir + "/" +self.option("sample_name"))
        os.mkdir(self.output_dir + "/" +self.option("sample_name"))
        all_fasta = {}
        with open(self.option("fasta").prop["path"],"r") as v:
            data = v.read()
            for i in data.split(">"):
                if i.strip():
                    all_fasta[i.split("\n")[0].split("\t")[0].split(" ")[0]] = "".join(i.strip().split("\n")[1:])
        with open(self.phigaro_tsv,"r") as f, open(self.output_dir + "/" + self.option("sample_name") + "/" + "phigaro_stat.xls","w") as t:
            t.write("Prophage ID\tSample Name\tLocation\tStart\tEnd\tLen (bp)\tVOG No.\tTaxonomy\tGC Content (%)\n")
            data1 = f.readlines()
            ph_num = 1
            for i in data1[1:]:
                ph_id = "PH" + str(ph_num).zfill(2)
                ph_num += 1
                gc_content = self.get_gc(all_fasta[i.strip().split("\t")[0]][int(i.strip().split("\t")[1]):int(i.strip().split("\t")[2])])
                t.write(ph_id+"\t"+self.option("sample_name")+"\t"+i.strip().split("\t")[0]+"\t"+i.strip().split("\t")[1]+"\t"+
                        i.strip().split("\t")[2]+"\t"+str(abs(int(i.strip().split("\t")[2])-int(i.strip().split("\t")[1])))+"\t"+
                        str(len(i.strip().split("\t")[-1].split(",")))+"\t"+ i.strip().split("\t")[4] +"\t"+str(gc_content)+"\n")

        with open(self.output_dir + "/" + self.option("sample_name") + "/" + "phigaro_detail.xls","w") as tt, open(self.phigaro_gff3,"r") as ff:
            tt.write("Prophage ID\tVOG ID\tSample Name\tLocation\tStart\tEnd\tLen (bp)\tGC Content (%)\n")
            data2 = ff.readlines()
            for line in data2:
                if line.startswith("#"):
                    pass
                else:
                    if line.strip().split("\t")[2] == "gene":
                        gc_content = self.get_gc(
                            all_fasta[line.split("\t")[0]][int(line.split("\t")[3]):int(line.split("\t")[4])])
                        prophage = "PH" + line.strip().split("\t")[-1].split(";")[1].split("prophage")[-1].zfill(2)
                        if line.strip().split("\t")[-1].split(";")[-1].split("=")[-1] == ".":
                            vog = "-"
                        else:
                            vog = line.strip().split("\t")[-1].split(";")[-1].split("=")[-1]
                            tt.write(
                                prophage + "\t" + vog + "\t" + self.option("sample_name") + "\t" + line.split("\t")[
                                    0] + "\t" +
                                line.split("\t")[3] + "\t" + line.split("\t")[4] + "\t" + str(
                                    abs(int(line.strip().split("\t")[4]) - int(line.strip().split("\t")[3])))
                                + "\t" + str(gc_content) + "\n")

    def set_output(self):
        """
        设置结果文件目录
        """
        if os.path.exists(self.output_dir + "/" + self.option("sample_name") + "/" + "phigaro.gff3"):
            os.remove(self.output_dir + "/" + self.option("sample_name") + "/" + "phigaro.gff3")
        os.link(self.phigaro_gff3,self.output_dir + "/" + self.option("sample_name") + "/" + "phigaro.gff3")
        self.logger.info("生成结果文件完成")

    def get_gc(self,str):
        """
        提取字符串的gc含量
        :param str:
        :return:
        """
        GC_len = (str.count('g') + str.count('G') + str.count('c') + str.count('C'))
        N_len = (str.count('n') + str.count('N'))
        AT_len = (str.count('a') + str.count('A') + str.count('t') + str.count('T'))
        all_len = GC_len + N_len + AT_len
        gc_percent = round(GC_len / float(all_len) * 100,2)
        return gc_percent

    def run(self):
        """
        运行
        """
        super(PhigaroTool, self).run()
        num = self.run_phigaro()
        if num<=1:
            self.end()
        else:
            self.run_stat()
            self.set_output()
            self.end()
# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
# version 1.0
# last_modify: 2021.04.07

import os, time
import re, shutil
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from mbio.packages.bac_comp_genome.common_function import link_dir
from Bio import SeqIO

class SrnaPredictAgent(Agent):
    """
    基因组sRNA预测小工具
    采用conda进行
    """
    def __init__(self, parent):
        super(SrnaPredictAgent, self).__init__(parent)
        options = [
            {"name": "fasta", "type": "infile", "format": "sequence.fasta"},# 序列文件
            {"name": "sample_name", "type": "string"},# 样本名
            {"name": "species", "type": "string", "default": "Bacteria"},# 物种 Bacteria,Archaea,Viruses
        ]
        self.add_option(options)
        #self._memory_increase_step = 20

    def check_options(self):
        if not self.option("fasta").is_set:
            raise OptionError("必须输入序列文件")

    def set_resource(self):
        self._cpu = 8
        self._memory = '40G'

    def end(self):
        super(SrnaPredictAgent, self).end()

class SrnaPredictTool(Tool):
    def __init__(self, config):
        super(SrnaPredictTool, self).__init__(config)
        self.seqstat = self.config.SOFTWARE_DIR + "/bioinfo/seq/biosquid_1.9g+cvs20050121/bin/seqstat"
        self.conda = self.config.SOFTWARE_DIR + "/program/miniconda3/"
        self.python_script = self.config.PACKAGE_DIR + "/bacgenome/"
        self.python = "/program/Python/bin/python"
        self.set_environ(PATH=self.conda + "bin")

    def seq_stat(self):
        """
        统计总长度
        :return:
        """
        all_length = 0
        for i in SeqIO.parse(self.option("fasta").prop["path"], "fasta"):
            all_length += len(i.seq)
        length = round(float(all_length) * 2 / 1000000, 5)
        self.logger.info(all_length)
        self.logger.info(length)
        return all_length,length

    def run_cmscan(self):
        """
        运行cmscan
        :return:
        """
        self.total,seq_length = self.seq_stat()
        self.logger.info(seq_length)
        cmd = "{} -Z {} -o {}_assembly.txt --tblout {}_assembly.tblout --fmt 2 --noali --cut_ga --rfam " \
              "--nohmmonly --cpu 8 {} {}".format("program/miniconda3/bin/cmscan", str(seq_length), self.option("sample_name"), self.option("sample_name"), self.conda + "db/cm/" + self.option("species"), self.option("fasta").prop["path"])
        self.logger.info(cmd)
        command = self.add_command("run_cmscan", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_cmscan运行完成！")
        else:
            self.set_error("run_cmscan运行出错!")

    def run_stat(self):
        """
        运行run_stat,主要进行数据统计，改名字，整理数据，提取序列
        :return:
        """
        nu = self.get_num(self.work_dir+"/"+ self.option("sample_name") + "_assembly.tblout")
        if nu >=1:
            cmd = "{} {}srna_stat.py --fa {} --in {} --s {} --len {} --o {}".format(self.python, self.python_script,
                                                                                    self.option("fasta").prop["path"],
                                                                                    self.work_dir + "/" + self.option(
                                                                                        "sample_name") + "_assembly.tblout",
                                                                                    self.option("sample_name"),
                                                                                    str(self.total),
                                                                                    self.output_dir + "/" + self.option(
                                                                                        "sample_name"))
            self.logger.info(cmd)
            command = self.add_command("run_stat", cmd).run()
            self.wait(command)
            if command.return_code == 0:
                self.logger.info("run_stat运行完成！")
            else:
                self.set_error("run_stat运行出错!")

    def get_num(self,input):
        with open(input, "r") as f:
            lines = f.readlines()
            num = 0
            for line in lines:
                if re.search("#", line):
                    continue
                else:
                    num +=1
            return num

    def run(self):
        """
        运行
        """
        super(SrnaPredictTool, self).run()
        self.run_cmscan()
        self.run_stat()
        self.end()
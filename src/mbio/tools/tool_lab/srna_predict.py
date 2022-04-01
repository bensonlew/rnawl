# -*- coding: utf-8 -*-
# __author__ = 'zzg'
# version 1.0
# last_modify: 2020.02.23

import os, time
import re, shutil
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from mbio.packages.bac_comp_genome.common_function import link_dir

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
        self.set_environ(PATH = self.conda +"bin")

    def seq_stat(self):
        """
        统计总长度
        :return:
        """
        length = 5
        tmp_data = os.popen(self.seqstat + " " + self.option("fasta").prop["path"])
        for i in tmp_data.read().split("\n"):
            if i.startswith("Total"):
                all_length = i.split(" ")[-1].strip()
                length = round(float(all_length) * 2 / 1000000,5)
                with open(self.work_dir + "/length_stat.txt", "w") as t:
                    t.write(all_length)
                break
        return length

    def run_cmscan(self):
        """
        运行cmscan
        :return:
        """
        seq_length = self.seq_stat()
        cmd = "{} -Z {} -o {}_assembly.txt --tblout {}_assembly.tblout --fmt 2 --noali --cut_ga --rfam " \
              "--nohmmonly --cpu 8 {} {}".format("program/miniconda3/bin/cmscan", str(seq_length), self.option("sample_name"), self.option("sample_name"), self.conda + "db/cm/" + self.option("species"), self.option("fasta").prop["path"])
        self.logger.info(cmd)
        command = self.add_command("run_cmscan", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_roary运行完成！")
        else:
            self.set_error("run_roary运行完成运行出错!")

    def set_output(self):
        """
        设置结果文件目录
        """
        outfile1 = os.path.join(self.output_dir, self.option("sample_name") + "_assembly.txt")
        outfile2 = os.path.join(self.output_dir, self.option("sample_name") + "_assembly.tblout")
        outfile3 = os.path.join(self.output_dir, "length_stat.txt")
        if os.path.exists(outfile1):
            os.remove(outfile1)
        if os.path.exists(outfile2):
            os.remove(outfile2)
        if os.path.exists(outfile3):
            os.remove(outfile3)
        os.link(self.work_dir + "/" + self.option("sample_name") + "_assembly.txt", outfile1)
        os.link(self.work_dir + "/" + self.option("sample_name") + "_assembly.tblout", outfile2)
        os.link(self.work_dir + "/length_stat.txt", outfile3)
        self.logger.info("生成结果文件完成")

    def run(self):
        """
        运行
        """
        super(SrnaPredictTool, self).run()
        #self.seq_stat()
        self.run_cmscan()
        self.set_output()
        self.end()
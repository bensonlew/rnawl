# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# last_modified: 20171130

"""fastq序列去掉前三个碱基，保留左边前51bp"""
import os
import glob
import re
from biocluster.config import Config
from biocluster.tool import Tool
from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError


class CutLeft51Agent(Agent):
    def __init__(self, parent):
        super(CutLeft51Agent, self).__init__(parent)
        options = [
            {"name": "fastq", "type": "infile", "format": "datasplit.fastq"},  # 样本fastq文件
            {"name": "fastq_l", "type": "infile", "format": "datasplit.fastq"},  # 样本左端fastq文件
            {"name": "fastq_r", "type": "infile", "format": "datasplit.fastq"},  # 样本右端fastq文件
            {"name": "cut_three", "type": "string", "default": "True"},  # 是否删除前三个碱基
            {"name": "cut_fastq", "type": "outfile", "format": "datasplit.fastq"},  # 去掉前3bp，取前51bp的fastq
            {"name": "cut_fastq_l", "type": "outfile", "format": "datasplit.fastq"},  # 去掉前3bp，取前51bp的fastq_l
            {"name": "cut_fastq_r", "type": "outfile", "format": "datasplit.fastq"}  # 去掉前3bp，取前51bp的fastq_r
        ]
        self.add_option(options)
        # self.queue = "chaifen"  # 投递到指定的队列chaifen

    def check_options(self):
        if self.option("fastq").is_set:
            pass
        elif self.option("fastq_l").is_set and self.option("fastq_r").is_set:
            pass
        else:
            raise OptionError("请设置fastq文件或者fastq_l和fastq_r文件")

    def set_resource(self):
        self._cpu = "2"
        self._memory = "10G"

    def end(self):
        super(CutLeft51Agent, self).end()


class CutLeft51Tool(Tool):
    def __init__(self, config):
        super(CutLeft51Tool, self).__init__(config)
        self._version = 1.0
        self.python = "program/Python/bin/python"
        self.cut_left = self.config.PACKAGE_DIR + "/datasplit/cut_left51_v1.py"

    def run_single_cut_left51(self):
        """单端序列删除序列的前三个碱基后保留前51bp"""
        fastq_file = self.option("fastq").prop["path"]
        outfile = os.path.join(self.output_dir, os.path.basename(fastq_file).split(".")[0] + ".fq")
        cmd = "{} {} -i {} -o {} -cut_three {}".format(self.python, self.cut_left, fastq_file, outfile, self.option("cut_three"))
        command = self.add_command("cut_left51", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("cut_lest51运行成功")
        else:
            self.set_error("cut_lest51运行失败")
        self.option("cut_fastq", outfile)

    def run_pariend_cut_left51(self):
        """双端序列删除序列的前三个碱基后保留前51bp"""
        fastq_file = self.option("fastq_l").prop["path"]
        outfile_l = os.path.join(self.output_dir, os.path.basename(fastq_file).split(".")[0] + ".fq")
        cmd = "{} {} -i {} -o {} -cut_three {}".format(self.python, self.cut_left, fastq_file, outfile_l, self.option("cut_three"))
        command = self.add_command("cut_left51_l", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("cut_lest51左端运行成功")
        else:
            self.set_error("cut_lest51左端运行失败")
        fastq_file = self.option("fastq_r").prop["path"]
        outfile_r = os.path.join(self.output_dir, os.path.basename(fastq_file).split(".")[0] + ".fq")
        cmd = "{} {} -i {} -o {} -cut_three {}".format(self.python, self.cut_left, fastq_file, outfile_r, self.option("cut_three"))
        command = self.add_command("cut_left51_r", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("cut_lest51右端运行成功")
        else:
            self.set_error("cut_lest51右端运行失败")
        self.option("cut_fastq_l", outfile_l)
        self.option("cut_fastq_r", outfile_r)

    def run(self):
        super(CutLeft51Tool, self).run()
        if self.option("fastq").is_set:
            self.run_single_cut_left51()
        else:
            self.run_pariend_cut_left51()
        self.end()

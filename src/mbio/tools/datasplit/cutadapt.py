# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# last_modified: 20171129

"""cutadapt 用于ncRNA质控"""
import os
from biocluster.config import Config
from biocluster.tool import Tool
from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError


class CutadaptAgent(Agent):
    def __init__(self, parent=None):
        super(CutadaptAgent, self).__init__(parent)
        options = [
            {"name": "fastq", "type": "infile", "format": "sequence.fastq"},  # 样本fastq文件
            {"name": "fastq_l", "type": "infile", "format": "sequence.fastq"},  # 样本左端fastq文件
            {"name": "fastq_r", "type": "infile", "format": "sequence.fastq"},  # 样本右端fastq文件
            {"name": "min_length", "type": "string", "default": "18"},  # reads最短长度设置，舍弃长度小于此值的序列
            {"name": "low_quality_base", "type": "string", "default": "33"},  # 过滤掉序列两端的低质量碱基数
            {"name": "adaptor", "type": "string", "default": "TGGAATTCTCGGGTGCCAAGG"},  # 单端接头
            {"name": "r_adaptor", "type": "string", "default": "TGGAATTCTCGGGTGCCAAGG"},  # 右端接头
            {"name": "l_adaptor", "type": "string", "default": "GATCGTCGGACTGTAGAACTCTGAAC"},  # 左端接头
        ]
        self.add_option(options)
        self.step.add_steps('cutadapt')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.cutadapt.start()
        self.step.update()

    def step_end(self):
        self.step.cutadapt.finish()
        self.step.update()

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


class CutadaptTool(Tool):
    def __init__(self, config):
        super(CutadaptTool, self).__init__(config)
        self._version = 1.0
        self.cutadapt = "program/Python/bin/cutadapt"

    def run_pairend_cutadapt(self):
        """双端fastq序列进行质控，cutadapt"""
        read1 = self.option("fastq_l").prop["path"]
        read2 = self.option("fastq_r").prop["path"]
        trim1 = os.path.join(self.output_dir, "sample_R1_clean.fastq")
        trim2 = os.path.join(self.output_dir,  "sample_R2_clean.fastq")
        cmd = "{} -q {} --minimum-length {} -a {} -A {} -o {} -p {} {} {}".format(self.cutadapt, \
               self.option("low_quality_base"), self.option("min_length"), self.option("r_adaptor"),\
               self.option("l_adaptor"), trim1, trim2, read1, read2)
        command = self.add_command("cutadapt", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("cutadapt双端完成")
        else:
            self.set_error("运行cutadapt双端出错")

    def run_single_cutadapt(self):
        """一个fastq进行cutadapt, 质控"""
        self.logger.info("开始运行cutadapt")
        read1 = self.option("fastq").prop["path"]
        trim1 = os.path.join(self.output_dir, "sample_clean.fastq")
        cmd = "{} -q {} --minimum-length {} -a {} -o {} {}".format(self.cutadapt, self.option("low_quality_base"),\
              self.option("min_length"), self.option("adaptor"), trim1, read1)
        command = self.add_command("cutadapt", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("运行cutadapt单端完成")
        else:
            self.set_error("运行cutadapt单端出错")

    def run(self):
        super(CutadaptTool, self).run()
        if self.option("fastq").is_set:
            self.run_single_cutadapt()
        else:
            self.run_pairend_cutadapt()
        self.end()

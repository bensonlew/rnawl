#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import division
from biocluster.agent import Agent
from biocluster.tool import Tool
import os
from biocluster.core.exceptions import OptionError
import glob
import re
import subprocess


class SeqPrepAgent(Agent):
    """
    seqprep:用于对PE序列去接头的工具
    version 1.0
    author: qindanhua
    """
    def __init__(self, parent):
        super(SeqPrepAgent, self).__init__(parent)
        options = [
            {"name": "fastq_r", "type": "infile", "format": "datasplit.fastq"},  # 输入文件PE的右端序列
            {"name": "fastq_l", "type": "infile", "format": "datasplit.fastq"},  # PE的左端序列
            {"name": "seqprep_r", "type": "outfile", "format": "datasplit.fastq"},  # PE的右端输出结果
            {"name": "seqprep_l", "type": "outfile", "format": "datasplit.fastq"},  # PE的左端输出结果
            {"name": "quality", "type": "int", "default": 30},
            {"name": "length", "type": "int", "default": 30},
            {"name": "adapter_a", "type": "string", "default": "AGATCGGAAGAGCACACGTC"},
            {"name": "adapter_b", "type": "string", "default": "AGATCGGAAGAGCGTCGTGT"},
        ]
        self.add_option(options)
        # self.queue = "chaifen"  # 投递到指定的队列chaifen

    def check_options(self):
        """
        检测参数是否正确
        """
        if not self.option("fastq_l").is_set:
            raise OptionError("缺少PE序列左端文件", code="31802402")
        if not self.option("fastq_r").is_set:
            raise OptionError("缺少PE序列右端文件", code="31802403")

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 1
        self._memory = ''

    def end(self):
        super(SeqPrepAgent, self).end()


class SeqPrepTool(Tool):
    """
    version 1.0
    """
    def __init__(self, config):
        super(SeqPrepTool, self).__init__(config)
        self.seqprep_path = '/bioinfo/seq/SeqPrep'
        if self.option("fastq_r").is_set and self.option("fastq_l").is_set:
            self.sample_name_r = self.option("fastq_r").prop["path"].split("/")[-1]
            self.sample_name_l = self.option("fastq_l").prop["path"].split("/")[-1]

    def seqprep(self):
        fq_r_path = self.option("fastq_r").prop['path']
        fq_l_path = self.option("fastq_l").prop['path']
        cmd = self.seqprep_path + " -f {} -r {} -1 {} -2 {} -q {} -L {} -A {} -B {}".\
            format(fq_l_path, fq_r_path, '{}_seqprep_l.gz'.format(self.sample_name_l), '{}_seqprep_r.gz'.format(self.sample_name_r), self.option('quality'),
                   self.option('length'), self.option("adapter_a"), self.option("adapter_b"))
        self.logger.info("开始运行seqprep")
        command = self.add_command("seqprep", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("运行seqprep完成")
        else:
            self.set_error("运行seqprep运行出错!")

    def set_output(self):
        """
        将结果文件链接至output
        """
        for f in os.listdir(self.output_dir):
            os.remove(os.path.join(self.output_dir, f))
        file_path = glob.glob(r"*.gz")
        for f in file_path:
            output_dir = os.path.join(self.output_dir, f)
            if os.path.exists(output_dir):
                os.remove(output_dir)
            os.link(os.path.join(self.work_dir, f), output_dir)
            os.remove(os.path.join(self.work_dir, f))
            if "seqprep_l.gz" in f:
                self.option("seqprep_l").set_path(os.path.join(self.output_dir, f))
            elif "seqprep_r.gz" in f:
                self.option("seqprep_r").set_path(os.path.join(self.output_dir, f))

    def run(self):
        """
        运行
        """
        super(SeqPrepTool, self).run()
        self.seqprep()
        self.set_output()
        self.end()

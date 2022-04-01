# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
# version 1.0
# last_modify: 2020.07.23

import os
import re, shutil
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError


class SplitFilesAgent(Agent):
    """
    主要是根据输入的宏基因组文件进行拆分处理
    """
    def __init__(self, parent):
        super(SplitFilesAgent, self).__init__(parent)
        options = [
            {"name": "genome_fa", "type": "infile", "format": "sequence.fasta"},  # 细菌的基因组fasta文件或宏基因组组装的fasta文件
            {"name": "gene_faa", "type": "infile",  "format": "sequence.fasta"},  # 预测的基因核酸文件
            {"name": "gene_fna", "type": "infile", "format": "sequence.fasta"},   # 预测的基因蛋白文件
            {"name": "gene_rnt", "type": "string"},  # rNA的统计文件
            {"name": "gene_gff", "type": "string"},   # 预测的基因gff文件
            {"name": "sample", "type": "string"}
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("genome_fa").is_set:
            raise OptionError("必须设置参数genome_fa文件!")
        if not self.option("gene_faa").is_set:
            raise OptionError("必须设置参数gene_faa文件!")
        if not self.option("gene_fna").is_set:
            raise OptionError("必须设置参数gene_fna文件!")
        if not self.option("gene_gff"):
            raise OptionError("必须设置参数gene_gff文件!")
        if not self.option("gene_rnt"):
            raise OptionError("必须设置参数gene_rnt文件!")

    def set_resource(self):
        size = os.path.getsize(self.option("genome_fa").prop['path'])
        num = float(size)/(1024*1024*1024)
        if num <5:
            self._cpu = 2
            self._memory = '10G'
        else:
            self._cpu = 2
            self._memory = str(int(num)*3) + "G"

    def end(self):
        super(SplitFilesAgent, self).end()

class SplitFilesTool(Tool):
    def __init__(self, config):
        super(SplitFilesTool, self).__init__(config)
        self.python = "/program/Python/bin/python"
        self.python_script = self.config.PACKAGE_DIR + "/mobile_genetic_elements/split_file.py"
        self.genome = self.option("genome_fa").prop['path']
        self.gene_fna = self.option("gene_fna").prop['path']
        self.gene_faa = self.option("gene_faa").prop['path']
        self.gene_gff = self.option("gene_gff")
        self.gene_rnt = self.option("gene_rnt")
        self.sample = self.option("sample")

    def run_file(self):
        cmd = "{} {} --s {} --a {} --f {} --n {} --g {} --t {} --o {}".format(self.python, self.python_script, self.genome, self.sample, self.gene_faa, self.gene_fna, self.gene_gff, self.gene_rnt, self.output_dir)
        command = self.add_command("run_file", cmd).run()
        self.wait(command)
        self.logger.info(cmd)
        if command.return_code == 0:
            self.logger.info("run_file运行完成！")
        else:
            self.set_error("run_file运行完成运行出错!")

    def run(self):
        super(SplitFilesTool, self).run()
        self.run_file()
        self.end()
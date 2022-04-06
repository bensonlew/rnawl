# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'
# version 1.0

import os
import re, shutil
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError


class IsGeneDirAgent(Agent):
    """
    主要是根据输入文件，处理成isescan程序所需要的输入文件
    """

    def __init__(self, parent):
        super(IsGeneDirAgent, self).__init__(parent)
        options = [
            {"name": "genome_fa", "type": "infile", "format": "sequence.fasta"},  # 细菌的基因组fasta文件或宏基因组组装的fasta文件
            {"name": "gene_faa", "type": "infile",  "format": "sequence.fasta"},  # 预测的基因蛋白文件
            {"name": "gene_fna", "type": "infile", "format": "sequence.fasta"},   # 预测的基因核酸文件
            {"name": "gene_gff", "type": "string"},   # 预测的基因gff文件
            {"name": "sample", "type": "string"}, ## 传入的样本名称
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

    def set_resource(self):
        self._cpu = 2
        self._memory = '10G'

    def end(self):
        super(IsGeneDirAgent, self).end()

class IsGeneDirTool(Tool):
    def __init__(self, config):
        super(IsGeneDirTool, self).__init__(config)
        self.python = "/miniconda2/bin/python"
        self.python_script = self.config.PACKAGE_DIR + "/bacgenome/is_gene_dir.py"
        self.genome = self.option("genome_fa").prop['path']
        self.gene_fna =self.option("gene_fna").prop['path']
        self.gene_faa =self.option("gene_faa").prop['path']
        self.gene_gff =self.option("gene_gff")

    def run_isfile(self):
        """
        format dir文件
        :return:
        """
        cmd = "{} {} --s {} --f {} --n {} --g {} --o {} --a {}".format(self.python, self.python_script, self.genome, self.gene_faa, self.gene_fna, self.gene_gff, self.output_dir, self.option("sample"))
        self.logger.info(cmd)
        command = self.add_command("run_isfile", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_isfile运行完成！")
        else:
            self.set_error("run_isfile运行完成运行出错!")

    def run(self):
        super(IsGeneDirTool, self).run()
        self.run_isfile()
        self.end()

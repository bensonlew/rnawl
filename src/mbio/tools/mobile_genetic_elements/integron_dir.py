# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
# version 1.0
# last_modify: 2020.06.12

import os
import re, shutil
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError


class IntegronDirAgent(Agent):
    """
    主要是根据输入文件，处理成Integron_finder程序所需要的输入文件
    """

    def __init__(self, parent):
        super(IntegronDirAgent, self).__init__(parent)
        options = [
            {"name": "genome_fa", "type": "infile", "format": "sequence.fasta"},  # 细菌的基因组fasta文件或宏基因组组装的fasta文件
            {"name": "gene_faa", "type": "infile",  "format": "sequence.fasta"},  # 预测的基因核酸文件
            {"name": "gene_gff", "type": "string"},   # 预测的基因gff文件
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("genome_fa").is_set:
            raise OptionError("必须设置参数genome_fa文件!")
        if not self.option("gene_faa").is_set:
            raise OptionError("必须设置参数gene_faa文件!")
        if not self.option("gene_gff"):
            raise OptionError("必须设置参数gene_gff文件!")

    def set_resource(self):
        self._cpu = 2
        self._memory = '10G'

    def end(self):
        super(IntegronDirAgent, self).end()

class IntegronDirTool(Tool):
    def __init__(self, config):
        super(IntegronDirTool, self).__init__(config)
        self.python = "/program/Python/bin/python"
        self.python_script = self.config.PACKAGE_DIR + "/mobile_genetic_elements/integron_dir.py"
        self.genome = self.option("genome_fa").prop['path']
        self.gene_faa =self.option("gene_faa").prop['path']
        self.gene_gff =self.option("gene_gff")

    def run_integron(self):
        cmd = "{} {} --s {} --f {} --g {} --o {}".format(self.python, self.python_script, self.genome, self.gene_faa, self.gene_gff, self.output_dir)
        command = self.add_command("run_isfile", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_integron运行完成！")
        else:
            self.set_error("run_integron运行完成运行出错!")

    def run(self):
        super(IntegronDirTool, self).run()
        self.run_integron()
        self.end()
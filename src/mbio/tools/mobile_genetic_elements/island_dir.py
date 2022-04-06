# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
# version 1.0
# last_modify: 2020.06.12

import os
import re, shutil
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError


class IslandDirAgent(Agent):
    """
    主要是根据输入文件，处理成Island_dimob程序所需要的输入文件
    """

    def __init__(self, parent):
        super(IslandDirAgent, self).__init__(parent)
        options = [
            {"name": "gene_fa", "type": "infile", "format": "sequence.fasta"},  # 细菌的基因核酸文件
            {"name": "gene_faa", "type": "infile",  "format": "sequence.fasta"},  # 预测的基因蛋白文件
            {"name": "gene_gff", "type": "string",},   # 预测的基因gff文件
            {"name": "genome_fa", "type": "infile", "format": "sequence.fasta"},
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("gene_fa").is_set:
            raise OptionError("必须设置参数gene_fa文件!")
        if not self.option("gene_faa").is_set:
            raise OptionError("必须设置参数gene_faa文件!")
        if not self.option("gene_gff"):
            raise OptionError("必须设置参数gene_gff文件!")

    def set_resource(self):
        self._cpu = 2
        self._memory = '10G'

    def end(self):
        super(IslandDirAgent, self).end()

class IslandDirTool(Tool):
    def __init__(self, config):
        super(IslandDirTool, self).__init__(config)
        self.python = "/miniconda2/bin/python"
        self.python_script = self.config.PACKAGE_DIR + "/mobile_genetic_elements/island_dir.py"
        self.genome = self.option("genome_fa").prop['path']
        self.gene = self.option("gene_fa").prop['path']
        self.gene_faa =self.option("gene_faa").prop['path']
        self.gene_gff =self.option("gene_gff")

    def run_island(self):
        cmd = "{} {} --s {} --f {} --n {} --g {} --o {}".format(self.python, self.python_script, self.genome, self.gene_faa, self.gene, self.gene_gff, self.output_dir)
        command = self.add_command("run_islandfile", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_island运行完成！")
        else:
            self.set_error("run_island运行完成运行出错!")

    def run(self):
        super(IslandDirTool, self).run()
        self.run_island()
        self.end()
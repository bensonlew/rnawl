# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
# version 1.0
# last_modify: 2020.06.12

import os
import re, shutil
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError


class GihunterDirAgent(Agent):
    """
    主要是根据输入文件，处理成Gihunter程序所需要的输入文件
    """

    def __init__(self, parent):
        super(GihunterDirAgent, self).__init__(parent)
        options = [
            {"name": "genome_fa", "type": "infile", "format": "sequence.fasta"},  # 细菌的基因组fasta文件或宏基因组组装的fasta文件
            {"name": "rnt", "type": "string"},  # 预测的基因核酸文件
            {"name": "gff", "type": "string"},   # 预测的基因gff文件
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("genome_fa").is_set:
            raise OptionError("必须设置参数genome_fa文件!")
        if not self.option("rnt"):
            raise OptionError("必须设置参数rnt文件!")
        if not self.option("gff"):
            raise OptionError("必须设置参数gff文件!")

    def set_resource(self):
        self._cpu = 2
        self._memory = '10G'

    def end(self):
        super(GihunterDirAgent, self).end()

class GihunterDirTool(Tool):
    def __init__(self, config):
        super(GihunterDirTool, self).__init__(config)
        self.python = "/miniconda2/bin/python"
        self.python_script = self.config.PACKAGE_DIR + "/mobile_genetic_elements/gihunter_dir.py"
        self.genome = self.option("genome_fa").prop['path']
        self.rnt =self.option("rnt")
        self.gff =self.option("gff")

    def run_gihunter(self):
        cmd = "{} {} --s {} --r {} --g {} --o {}".format(self.python, self.python_script, self.genome, self.rnt, self.gff, self.output_dir)
        command = self.add_command("run_gihunter", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_gihunter运行完成！")
        else:
            self.set_error("run_gihunter运行完成运行出错!")

    def run(self):
        super(GihunterDirTool, self).run()
        self.run_gihunter()
        self.end()
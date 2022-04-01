# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
# version 1.0
# last_modify: 20119.05.09

import os
import re, shutil
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from mbio.packages.bac_com_genome.common_function import link_dir

class PanInputAgent(Agent):
    """
    细菌比较基因组PGAP
    """
    def __init__(self, parent):
        super(PanInputAgent, self).__init__(parent)
        options = [
            {"name": "strains", "type": "sting"},  # 物种+物种
            {"name": "infile_dir", "type": "infile", "format": "bac_comp_genome.input_dir"},  # 输入文件夹
            {"name": "out", "type": "outfile", "format": "bac_comp_genome.pagp_dir"},#输出PAGP所需要格式数据
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("strains"):
            raise OptionError("必须设置参数strains的!")
        if not self.option("infile_dir").is_set:
            raise OptionError("必须设置参数infile_dir文件夹!")

    def set_resource(self):
        self._cpu = 2
        self._memory = '10G'

    def end(self):
        super(PanInputAgent, self).end()

class PanInputTool(Tool):
    def __init__(self, config):
        super(PanInputTool, self).__init__(config)
        self.strains = self.option("strains")
        self.fasta = self.option("infile_dir").prop['path']
        self.perl_path = "/program/perl/perls/perl-5.24.0/bin/perl"
        self.perl_script = self.config.SOFTWARE_DIR + "/bioinfo/compare_genome/software/PGAP-1.2.1/Converter_finished.pl"
        self.out = self.work_dir + '/out_result'

    def run_input(self):
        cmd = "{} {} -S {} -I {} -O {}".format(self.perl_path, self.perl_script, self.strains, self.fasta, self.out)
        command = self.add_command("run_input", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_input运行完成！")
        else:
            self.set_error("run_input运行完成运行出错!")

    def set_output(self):
        link_dir(self.out, self.output_dir)
        self.option("out", self.output_dir)

    def run(self):
        super(PanInputTool, self).run()
        self.run_pagp()
        self.set_output()
        self.end()

#!/usr/bin/env python
# -*- coding: utf-8 -*-
# __author__ = 'gaohao'

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os,shutil

class DastoolInputAgent(Agent):
    """
    用于DAS_toools的输入文件生成
    version 1.0
    author: gaohao
    last_modify: 2019.01.08
    """

    def __init__(self, parent):
        super(DastoolInputAgent, self).__init__(parent)
        options = [
            {"name": "sof", "type": "string"},  #软件名称
            {"name": "bin_dir", "type": "infile", "format": "sequence.fasta_dir"},  # bin的目录
            {"name": "out", "type": "outfile", "format": "sequence.profile_table"},#生成bin的目录
        ]
        self.add_option(options)

    def check_options(self):
        """
        检测参数是否正确
        """
        if not self.option('bin_dir').is_set:
            raise OptionError("bind的文件目录不存在！")

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 2
        self._memory ="5G"

    def end(self):
        super(DastoolInputAgent, self).end()


class DastoolInputTool(Tool):
    """
    version 1.0
    """
    def __init__(self, config):
        super(DastoolInputTool, self).__init__(config)
        self.path =self.config.SOFTWARE_DIR + "/program/perl-5.24.0/bin"
        self.set_environ(PATH=self.path)
        self.Scaffolds2Bin_sh ="../../../../../.." + self.config.PACKAGE_DIR + '/metagbin/fasta_to_Scaffolds2Bin.sh'
        self.run_sh =self.config.SOFTWARE_DIR + "/bioinfo/metaGenomic/DAS_Tool-1.1.0/src/Fasta_to_Scaffolds2Bin.sh"
        self.dir =self.option('bin_dir').prop['path']
        self.sof = self.option('sof')
        self.out =self.work_dir + '/' + self.sof + ".scaffolds2bin.tsv"

    def run_input(self):
        cmd = "{} {} {} {} {}".format(self.Scaffolds2Bin_sh,self.run_sh,self.dir,self.sof,self.out)
        self.logger.info(cmd)
        self.logger.info("开始运行run_input")
        command = self.add_command("run_input", cmd)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("运行run_input完成")
        else:
            self.set_error("运行run_input运行出错!")

    def set_output(self):
        if os.path.exists(self.output_dir + '/' + self.sof + ".scaffolds2bin.tsv"):
            shutil.rmtree(self.output_dir + '/' + self.sof + ".scaffolds2bin.tsv")
        os.link(self.work_dir + '/' + self.sof + ".scaffolds2bin.tsv",self.output_dir + '/' + self.sof + ".scaffolds2bin.tsv")
        self.option('out',self.output_dir + '/' + self.sof + ".scaffolds2bin.tsv")

    def run(self):
        """
        运行
        """
        super(DastoolInputTool, self).run()
        self.run_input()
        self.set_output()
        self.end()
#!/usr/bin/env python
# -*- coding: utf-8 -*-
# __author__ = 'gaohao'

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os,shutil
from mbio.packages.metagbin.common_function import link_dir,bin_rename

class GetBinFastaAgent(Agent):
    """
    提取新bin的序列
    version 1.0
    author: gaohao
    last_modify: 2019.04.09
    """

    def __init__(self, parent):
        super(GetBinFastaAgent, self).__init__(parent)
        options = [
            {"name": "bin_path", "type": "infile", "format": "sequence.fasta"},  # bin的文件
            {"name": 'bin_list', "type": "string"},
            {"name": 'new_bin', "type": "string"},
            {"name": "out", "type": "outfile", "format": "sequence.fasta"},  # 生成bin的目录
        ]
        self.add_option(options)

    def check_options(self):
        """
        检测参数是否正确
        """
        if not self.option('bin_path').is_set:
            raise OptionError("bin序列文件不存在！")
        if self.option('new_bin') == "":
            raise OptionError("bin的新改名称不存在！")

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 2
        self._memory ="5G"

    def end(self):
        super(GetBinFastaAgent, self).end()

class GetBinFastaTool(Tool):
    """
    version 1.0
    """
    def __init__(self, config):
        super(GetBinFastaTool, self).__init__(config)
        self.perl ="/program/perl/perls/perl-5.24.0/bin/perl"
        self.perl_path = self.config.PACKAGE_DIR + "/metagbin/get_bin_fasta.pl"
        self.fa = self.option('bin_path').path
        self.prefix = self.option('new_bin')
        self.binlist = self.option('bin_list')

    def run_getfasta(self):
        cmd = "{} {} {} {} {}".format(self.perl,self.perl_path,self.fa,'"' + self.binlist + '"',self.prefix)
        self.logger.info(cmd)
        self.logger.info("开始运行run_getfasta")
        command = self.add_command("run_getfasta", cmd)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("运行run_getfasta完成")
        else:
            self.set_error("运行run_getfasta运行出错!")

    def set_output(self):
        if os.path.exists(self.output_dir + '/' + self.prefix + '.fa'):
            shutil.rmtree(self.output_dir + '/' + self.prefix + '.fa')
        os.link(self.work_dir + '/' + self.prefix + '.fa',self.output_dir + '/' + self.prefix + '.fa')
        self.option('out',self.output_dir + '/' + self.prefix + '.fa')

    def run(self):
        """
        运行
        """
        super(GetBinFastaTool, self).run()
        self.run_getfasta()
        self.set_output()
        self.end()
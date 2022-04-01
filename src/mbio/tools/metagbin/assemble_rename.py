# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
# version 1.0
# last_modify: 20119.01.14

import os
import re
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError

class AssembleRenameAgent(Agent):
    """
    组装文件改名字
    """
    def __init__(self, parent):
        super(AssembleRenameAgent, self).__init__(parent)
        options = [
            {"name": "assmle_fa", "type": "infile", "format": "sequence.fasta"},  # 组装文件
            {"name": "sample_name", "type": "string"},
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("assmle_fa").is_set:
            raise OptionError("必须设置参数assmle_fa序列文件!")

    def set_resource(self):
        self._cpu = 2
        self._memory = '5G'

    def end(self):
        super(AssembleRenameAgent, self).end()

class AssembleRenameTool(Tool):
    def __init__(self, config):
        super(AssembleRenameTool, self).__init__(config)
        self.fasta = self.option("assmle_fa").prop['path']
        self.perl_path = "/program/perl-5.24.0/bin/perl"
        self.rename =self.config.PACKAGE_DIR + "/metagbin/assemble_rename.pl"
        self.cut = self.config.PACKAGE_DIR + "/metagbin/tiqu_seq.pl"

    def run_rename(self):
        if os.path.exists(self.work_dir + '/' + self.option('sample_name') +'.scaf.fa'):
            os.remove(self.work_dir + '/' + self.option('sample_name') +'.scaf.fa')
        cmd = '{} {} {} {}'.format(self.perl_path, self.rename,self.fasta,self.option('sample_name'))
        command = self.add_command("run_rename", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("assemble组装文件改名运行完成!")
        else:
            self.set_error("assemble组装文件改名运行出错!")

    def split_seq(self):
        cmd ='{} {} {} {}'.format(self.perl_path,self.cut,self.work_dir + '/' + self.option('sample_name') +'.scaf.fa',self.work_dir + '/' + self.option('sample_name') +'.scaf_1000.fa')
        command = self.add_command("split_seq", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("split_seq去除1000以下运行完成!")
            self.set_output()
        else:
            self.set_error("split_seq去除1000以下运行出错!")

    def set_output(self):
        if os.path.getsize(self.work_dir + '/' + self.option('sample_name') +'.scaf_1000.fa') > 0:
            os.link(self.work_dir + '/' + self.option('sample_name') +'.scaf_1000.fa',self.output_dir + '/' + self.option('sample_name') +'.scaf.fa')

    def run(self):
        super(AssembleRenameTool, self).run()
        self.run_rename()
        self.split_seq()
        self.end()
# -*- coding: utf-8 -*-
# __author__ = 'hao.gao'

import os, sys
# import shutil
from biocluster.core.exceptions import OptionError
from biocluster.agent import Agent
from biocluster.tool import Tool


class ScafSelectAgent(Agent):
    """
    细菌基因组的scaffold的序列根据N50筛选
    version: v1
    author: hao.gao
    last_modify: 2018.01.29
    """

    def __init__(self, parent):
        super(ScafSelectAgent, self).__init__(parent)
        options = [
            {"name": "seq_dir", "type": "infile", "format": "sequence.fasta_dir"},  # 组装序列文件夹
            {"name": "scf_seq", "type": "outfile", "format": "sequence.fasta"}  # 输出最佳组装scf文件
        ]
        self.add_option(options)

    def check_options(self):
        """
        检查参数
        :return:
        """
        if not self.option('seq_dir'):
            raise OptionError('必须输入seq_dir文件', code="31301601")

    def set_resource(self):
        """
        :return:
        """
        self._cpu = 2
        self._memory = "5G"

    def end(self):
        super(ScafSelectAgent, self).end()

class ScafSelectTool(Tool):
    def __init__(self, config):
        super(ScafSelectTool, self).__init__(config)
        self.dir = self.option('seq_dir').prop['path']
        self.perl_path = "/program/perl-5.24.0/bin/perl"
        self.perl_script = self.config.PACKAGE_DIR + "/bacgenome/"
        self.out = self.work_dir + "/N50.last.result"

    def run(self):
        """
        运行
        :return:
        """
        super(ScafSelectTool, self).run()
        self.select_scaf()
        self.set_output()
        self.end()

    def select_scaf(self):
        cmd = '{} {}N50_select.pl {} {}'.format(self.perl_path, self.perl_script, self.dir,self.out)
        command = self.add_command('select_scaf', cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("select_scaf运行完成" )
        else:
            self.set_error("select_scaf运行出错!", code="31301601")

    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        if os.path.exists(self.out):
            with open(self.out, 'r') as f:
                lines = f.readlines()
                line = lines[0].strip('\r\n').split('\t')[0]
                self.option('scf_seq').set_path(line)
            self.logger.info("设置最佳N50组装结果文件")


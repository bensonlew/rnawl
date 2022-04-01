#!/usr/bin/env python
# -*- coding: utf-8 -*-
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os,re


class IslandIslanderAgent(Agent):
    """
    使用islander预测基因组岛
    version 1.0
    author: gaohao
    last_modify: 2020.07.02
    """

    def __init__(self, parent):
        super(IslandIslanderAgent, self).__init__(parent)
        options = [
            {"name": "fna", "type": "infile", "format": "sequence.fasta"},  #
            {"name": "out", "type": "outfile", "format": "bacgenome.island"},
        ]
        self.add_option(options)
        self.list =[]


    def check_options(self):
        """
        检测参数是否正确
        """
        if not self.option('fna').is_set:
            raise OptionError("请设置基因组序列文件！")

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 5
        self._memory = '20G'

    def end(self):
        super(IslandIslanderAgent, self).end()


class IslandIslanderTool(Tool):
    """
    version 1.0
    """
    def __init__(self, config):
        super(IslandIslanderTool, self).__init__(config)
        self.perl_path = "/program/perl-5.24.0/bin/perl"
        self.perl5path = self.config.SOFTWARE_DIR + "/bioinfo/Genomic/Sofware/tRNAscan-SE-2.0/lib/tRNAscan-SE:" + self.config.SOFTWARE_DIR + "/program/perl-5.24.0/lib"
        self.set_environ(PERL5LIB=self.perl5path)
        self.islander = self.config.SOFTWARE_DIR + "/bioinfo/Genomic/Sofware/Islander_software/Islander.pl"
        self.sample = os.path.basename(self.option("fna").prop['path'])

    def run_islander(self):
        if os.path.exists(self.work_dir + "/" + self.sample):
            os.remove(self.work_dir + "/" + self.sample)
        os.link(self.option("fna").prop['path'], self.work_dir + "/" + self.sample)
        fa = self.work_dir + "/" + self.sample.split('.fna')[0]
        cmd = "{} {} {} --verbose --translate --trna --annotate --reisland --table 11 --nocheck".format(self.perl_path, self.islander, fa)
        self.logger.info(cmd)
        path = 'run_islander'
        self.logger.info("开始运行%s " % path)
        command = self.add_command(path, cmd)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("运行%s完成" % path)
        else:
            self.set_error("运行%s运行出错!" % path)

    def set_output(self):
        path = self.work_dir + '/all_final_arrayed_islands.gff'
        if os.path.exists(self.work_dir + '/all_final_arrayed_islands.gff'):
            self.option('out').set_path(path)

    def run(self):
        """
        运行
        """
        super(IslandIslanderTool, self).run()
        self.run_islander()
        self.set_output()
        self.end()

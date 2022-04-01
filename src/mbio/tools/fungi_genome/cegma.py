# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
# version 1.0
# last_modify: 2018.05.28

import os
import re
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError

class CegmaAgent(Agent):
    """
    真菌基因组cemga基因组评估
    """
    def __init__(self, parent):
        super(CegmaAgent, self).__init__(parent)
        options = [
            {"name": "scaf_fa", "type": "infile", "format": "sequence.fasta"},  ##计算时使用base数量
            {"name": "sample_name", "type": "string"},
            ]
        self.add_option(options)

    def check_options(self):
        if not self.option("scaf_fa").is_set:
            raise OptionError("必须添加scaf_fa的fa文件！", code="32100501")

    def set_resource(self):
        self._cpu = 4
        self._memory = '50G'

    def end(self):
        super(CegmaAgent, self).end()

class CegmaTool(Tool):
    def __init__(self, config):
        super(CegmaTool, self).__init__(config)
        self.scaf = self.option("scaf_fa").prop['path']
        self.sample = self.option("sample_name")
        self.set_environ(WISECONFIGDIR=self.config.SOFTWARE_DIR + '/bioinfo/Genomic/Sofware/wise-2.4.1/wisecfg')
        self.set_environ(PATH=self.config.SOFTWARE_DIR + '/program/perl-5.24.0/bin:' + self.config.SOFTWARE_DIR + '/bioinfo/Genomic/Sofware/wise-2.4.1/src/bin:' + self.config.SOFTWARE_DIR + '/bioinfo/Genomic/Sofware/geneid/bin:' + self.config.SOFTWARE_DIR + '/bioinfo/Genomic/Sofware/ncbi-blast-2.2.28+/bin:' + self.config.SOFTWARE_DIR + '/bioinfo/Genomic/Sofware/hmmer-3.0/bin:' + self.config.SOFTWARE_DIR + '/bioinfo/Genomic/Sofware/CEGMA_v2-master/bin')
        self.set_environ(CEGMA=self.config.SOFTWARE_DIR + '/bioinfo/Genomic/Sofware/CEGMA_v2-master')
        # for ningbo guanqing.zou 20181023
        self.set_environ(PERL5LIB=self.config.SOFTWARE_DIR + '/bioinfo/Genomic/Sofware/CEGMA_v2-master/lib:' + self.config.SOFTWARE_DIR + '/program/perl-5.24.0/lib/5.24.0/x86_64-linux-thread-multi')  # add
        self.perl_path = "/program/perl-5.24.0/bin/perl"
        self.perl_script = self.config.PACKAGE_DIR + "/fungi_genome/"
        self.cegma = "bioinfo/Genomic/Sofware/CEGMA_v2-master/bin/cegma"
        self.file = self.work_dir + '/' + self.sample + '.completeness_report'

    def run_cegma(self):
        cmd = '{} -g {} -o {}'.format(self.cegma,self.scaf,self.sample)
        command = self.add_command("run_cegma", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_cegma运行完成")
        else:
            self.set_error("run_cegma运行出错!", code="32100501")

    def run_stat(self):
        cmd = '{} {}tiqu_cegma.pl {} {}'.format(self.perl_path,self.perl_script,self.file,self.work_dir + '/' + self.sample + '_cegma.xls')
        command = self.add_command("run_stat", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_stat运行完成")
        else:
            self.set_error("run_stat运行出错!", code="32100502")

    def set_output(self):
        pass

    def run(self):
        super(CegmaTool, self).run()
        self.run_cegma()
        self.run_stat()
        self.set_output()
        self.end()
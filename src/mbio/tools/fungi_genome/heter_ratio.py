# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
# version 1.0
# last_modify: 2018.05.25

import os
import re
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError

class HeterRatioAgent(Agent):
    """
    真菌基因组杂合率评估
    """
    def __init__(self, parent):
        super(HeterRatioAgent, self).__init__(parent)
        options = [
            {"name": "scaf_fa", "type": "infile", "format": "sequence.fasta"},  ##计算时使用base数量
            {"name": "heter", "type": "string"},  ##计算时使用base数量
            ]
        self.add_option(options)

    def check_options(self):
        if not self.option("scaf_fa").is_set:
            raise OptionError("必须添加scaf_fa的fa文件！", code="32101201")
        if not self.option("heter"):
            raise OptionError("必须添加heter的参数！", code="32101202")

    def set_resource(self):
        self._cpu = 2
        self._memory = '10G'

    def end(self):
        super(HeterRatioAgent, self).end()

class HeterRatioTool(Tool):
    def __init__(self, config):
        super(HeterRatioTool, self).__init__(config)
        self.scaf = self.option("scaf_fa").prop['path']
        self.heter = self.option("heter")
        self.sh_path = "../../../../../.." + self.config.PACKAGE_DIR + '/sequence/scripts/'
        self.set_environ(LD_LIBRARY_PATH=self.config.SOFTWARE_DIR + '/gcc/5.1.0/lib64')
        self.perl_path = "/program/perl-5.24.0/bin/perl"
        self.perl_drop = self.config.PACKAGE_DIR + "/sequence/scripts/"
        self.readsim = "bioinfo/Genomic/Sofware/simulate_heterozygosis/readsim"
        self.kmerfreq =self.config.SOFTWARE_DIR + "/bioinfo/Genomic/Sofware/kmerfreq_v5.0_BigKmer/Big_kmerfreq"
        self.read_list = self.work_dir + "/" + 'read.list'
        self.kmerfreq_file = self.work_dir + "/" + self.heter + '_kmer.freq'
        self.drop_file = self.work_dir + '/all.dropN.fasta'

    def run_drop_n(self):
        cmd = '{} {}drop_n_seq.pl {} {}'.format(self.perl_path,self.perl_drop,self.scaf,self.drop_file)
        command = self.add_command("run_drop_n", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_drop_n运行完成")
        else:
            self.set_error("run_drop_n运行出错!", code="32101201")

    def run_readsim(self):
        cmd = '{} {} -R 125 -i 500 -d 25 -X 30 -f 0.001 -H {}'.format(self.readsim,self.drop_file,self.heter)
        command = self.add_command("run_readsim", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_readsim运行完成")
        else:
            self.set_error("run_readsim运行出错!", code="32101202")

    def run_kmerfreq(self):
        os.system('ls read_*fq > %s' % self.read_list)
        cmd = '{}kmerfreq.sh {} {} {}'.format(self.sh_path,self.kmerfreq,self.read_list,self.kmerfreq_file)
        command = self.add_command("run_kmerfreq", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_kmerfreq运行完成")
        else:
            self.set_error("run_kmerfreq运行出错!", code="32101203")

    def set_output(self):
        if os.path.exists(self.output_dir + "/" + self.heter + '_kmer.freq'):
            os.remove(self.output_dir + "/" + self.heter + '_kmer.freq')
        os.link(self.kmerfreq_file,self.output_dir + "/" + self.heter + '_kmer.freq')

    def run(self):
        super(HeterRatioTool, self).run()
        self.run_drop_n()
        self.run_readsim()
        self.run_kmerfreq()
        self.set_output()
        self.end()
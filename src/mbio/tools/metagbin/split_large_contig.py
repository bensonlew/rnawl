 # -*- coding: utf-8 -*-
# __author__ = 'gaohao'
# version 1.0
# last_modify: 20119.01.21

import os
import re
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError

class SplitLargeContigAgent(Agent):
    """
    对组装的数据进行拆分，大于655360bp的contigs不能进行cd-hit聚类
    """
    def __init__(self, parent):
        super(SplitLargeContigAgent, self).__init__(parent)
        options = [
            {"name": "sca_fa", "type": "infile", "format": "sequence.fasta"},  # 每个bin注释的*.anno.xls文件目录
            {"name":"out_large","type":"outfile","format":"sequence.fasta"},
            {"name": "out", "type": "outfile", "format": "sequence.fasta"},
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("sca_fa").is_set:
            raise OptionError("必须设置参数sca_fa的文件!")

    def set_resource(self):
        self._cpu = 2
        self._memory = '5G'

    def end(self):
        super(SplitLargeContigAgent, self).end()

class SplitLargeContigTool(Tool):
    def __init__(self, config):
        super(SplitLargeContigTool, self).__init__(config)
        self.fasta = self.option("sca_fa").prop['path']
        self.perl_path = "/program/perl-5.24.0/bin/perl"
        self.perl_script =self.config.PACKAGE_DIR + "/metagbin/split_large_contig.pl"
        self.result = 'all'

    def run_split(self):
        cmd = '{} {} {} {}'.format(self.perl_path, self.perl_script,self.fasta,self.result)
        command = self.add_command("run_split", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("scaffolds拆分运行完成!")
        else:
            self.set_error("scaffolds拆分运行出错!")

    def set_output(self):
        for i in ['all.large.scaf.fa','all.scaf.fa']:
            if os.path.exists(self.output_dir + '/' + i):
                os.remove(self.output_dir + '/' + i)
            os.link(self.work_dir + '/' + i, self.output_dir + '/' + i)
        if os.path.getsize(self.work_dir + '/all.large.scaf.fa') >0:
            self.option('out_large', self.output_dir + '/all.large.scaf.fa')
        self.option('out', self.output_dir + '/all.scaf.fa')

    def run(self):
        super(SplitLargeContigTool, self).run()
        self.run_split()
        self.set_output()
        self.end()
#!/usr/bin/env python
# -*- coding: utf-8 -*-
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os,re
from mbio.packages.bac_comp_genome.common_function import link_dir,link_file

class AntismashFastaAgent(Agent):
    """
    前噬菌体的序列提取
    version 1.0
    author: gaohao
    last_modify: 2019.10.10
    """

    def __init__(self, parent):
        super(AntismashFastaAgent, self).__init__(parent)
        options = [
            {"name": "sample_name", "type": "string"},
            {"name": "antismash", "type": "string"},
            {"name": "genome", "type": "infile", "format": "sequence.fasta"},
        ]
        self.add_option(options)


    def check_options(self):
        """
        检测参数是否正确
        """
        if not self.option('antismash'):
            raise OptionError("请设置antismash的结果文件！")
        if not self.option('genome').is_set:
            raise OptionError("请设置genome软件的结果文件！")
        if not self.option('sample_name'):
            raise OptionError("请设置样品名！")

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 2
        self._memory = '5G'

    def end(self):
        super(AntismashFastaAgent, self).end()


class AntismashFastaTool(Tool):
    """
    version 1.0
    """
    def __init__(self, config):
        super(AntismashFastaTool, self).__init__(config)
        self.antismash = self.option("antismash")
        self.scaf_seq = self.option("genome").prop["path"]
        self.sample = self.option('sample_name')
        self.perl_path = "/program/perl-5.24.0/bin/perl"
        self.perl_script = self.config.PACKAGE_DIR + "/bac_comp_genome/"

    def tiqu_fasta(self):
        if os.path.exists(self.work_dir + '/' + self.sample + '_antismash.fna'):
            os.remove(self.work_dir + '/' + self.sample + '_antismash.fna')
        cmd = '{} {}get_antismash_fasta.pl {} {} {}'.format(self.perl_path, self.perl_script, self.scaf_seq, self.antismash, self.work_dir + '/' + self.sample + '_antismash.fna')
        self.logger.info(cmd)
        command = self.add_command("tiqu_fasta", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("tiqu_fasta运行完成!")
        else:
            self.set_error("tiqu_fasta运行出错!")

    def set_output(self):
        link_file(self.work_dir + '/' + self.sample + '_antismash.fna', self.output_dir + '/' + self.sample + '_antismash.fna')

    def run(self):
        """
        运行
        """
        super(AntismashFastaTool, self).run()
        self.tiqu_fasta()
        self.set_output()
        self.end()
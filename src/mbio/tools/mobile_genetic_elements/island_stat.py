#!/usr/bin/env python
# -*- coding: utf-8 -*-
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os,re
from mbio.packages.bac_comp_genome.common_function import link_dir,link_file

class IslandStatAgent(Agent):
    """
    生成基因组岛文件
    version 1.0
    author: gaohao
    last_modify: 2020.10.14
    """
    def __init__(self, parent):
        super(IslandStatAgent, self).__init__(parent)
        options = [
            {"name": "dir", "type": "string"},  #
            {"name": "sample_name", "type": "string"},
            {"name": "gff", "type": "string"},
            {"name": "genome", "type": "infile", "format": "sequence.fasta"},
        ]
        self.add_option(options)


    def check_options(self):
        """
        检测参数是否正确
        """
        if not self.option('dir'):
            raise OptionError("请设置dir文件目录！")
        if not self.option('genome').is_set:
            raise OptionError("请设置gemone序列文件！")
        if not self.option('sample_name'):
            raise OptionError("请设置样品名！")
        if not self.option('gff'):
            raise OptionError("请设置输入文件gff！")

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 2
        self._memory = '5G'

    def end(self):
        super(IslandStatAgent, self).end()


class IslandStatTool(Tool):
    """
    version 1.0
    """
    def __init__(self, config):
        super(IslandStatTool, self).__init__(config)
        self.gff = self.option("gff")
        self.scaf_seq = self.option("genome").prop["path"]
        self.dir =self.option('dir')
        self.sample = self.option('sample_name')
        self.python_path = "/miniconda2/bin/python"
        self.perl_path = "/program/perl-5.24.0/bin/perl"
        self.script = self.config.PACKAGE_DIR + "/mobile_genetic_elements/"

    def run_island(self):
        cmd = '{} {}combin_island.py --n {} --d {} --o {}'.format(self.python_path, self.script, self.sample, self.dir, self.work_dir + "/"+ self.sample +".island.xls")
        self.logger.info(cmd)
        command = self.add_command("run_island", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_island运行完成")
        else:
            self.set_error("run_island运行出错!")

    def run_island_stat(self):
        out = self.work_dir + "/"+ self.sample +".island.xls"
        cmd = '{} {}get_island.pl {} {} {}'.format(self.perl_path, self.script,
                                                           out, self.gff,self.sample)
        self.logger.info(cmd)
        command = self.add_command("run_island_stat", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_island_stat运行完成")
        else:
            self.set_error("run_island_stat运行出错!")

    def tiqu_fasta(self):
        if os.path.exists(self.work_dir + '/' + self.sample + '_island.fna'):
            os.remove(self.work_dir + '/' + self.sample + '_island.fna')
        cmd = '{} {}get_phage_fasta.pl {} {} {} {} {}'.format(self.perl_path, self.script, self.sample, "island",  self.scaf_seq,
                                                        self.work_dir + '/' + self.sample + '.GI_summary.xls',
                                                         self.work_dir + '/' + self.sample + '_island.fna')
        self.logger.info(cmd)
        command = self.add_command("tiqu_fasta", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("tiqu_fasta运行完成!")
        else:
            self.set_error("tiqu_fasta运行出错!")

    def set_output(self):
        link_file(self.work_dir + '/' + self.sample + '.GI_summary.xls', self.output_dir + '/' + self.sample + '.GI_summary.xls')
        link_file(self.work_dir + '/' + self.sample + '_island.fna', self.output_dir + '/' + self.sample + '_island.fna')

    def run(self):
        """
        运行
        """
        super(IslandStatTool, self).run()
        self.run_island()
        n = self.get_num(self.work_dir + "/"+ self.sample +".island.xls")
        if n >1:
            self.run_island_stat()
            self.tiqu_fasta()
            self.set_output()
            self.end()
        else:
            self.end()

    def get_num(self,file):
        with open(file,'r') as f:
            lines =f.readlines()
            num =len(lines)
        return num
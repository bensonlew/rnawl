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
    last_modify: 2018.04.13
    """

    def __init__(self, parent):
        super(IslandStatAgent, self).__init__(parent)
        options = [
            {"name": "diomb", "type": "infile", "format": "sequence.profile_table"},  #
            {"name": "islander", "type": "infile", "format": "bacgenome.island"},
            {"name": "sample_name", "type": "string"},
            {"name": "gff", "type": "infile", "format": "gene_structure.gff3"},
            {"name": "genome", "type": "infile", "format": "sequence.fasta"},
        ]
        self.add_option(options)


    def check_options(self):
        """
        检测参数是否正确
        """
        if not self.option('diomb').is_set:
            raise OptionError("请设置diomb软件的结果文件！")
        if not self.option('islander').is_set:
            raise OptionError("请设置islander软件的结果文件！")
        if not self.option('sample_name'):
            raise OptionError("请设置样品名！")
        if not self.option('gff').is_set:
            raise OptionError("请设置输入文件gff！")

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 2
        self._memory = '5G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            ["", "", ""],
            ["", "", ""]
        ])
        super(IslandStatAgent, self).end()


class IslandStatTool(Tool):
    """
    version 1.0
    """
    def __init__(self, config):
        super(IslandStatTool, self).__init__(config)
        self.gff = self.option("gff").prop["path"]
        self.scaf_seq = self.option("genome").prop["path"]
        self.diomb =self.option('diomb').prop['path']
        self.islander = self.option('islander').prop['path']
        self.sample = self.option('sample_name')
        self.perl_path = "/program/perl-5.24.0/bin/perl"
        self.perl_script = self.config.PACKAGE_DIR + "/bac_comp_genome/"

    def run_island(self):
        cmd = '{} {}combin_island.pl {} {}'.format(self.perl_path, self.perl_script, self.diomb,self.islander)
        self.logger.info(cmd)
        command = self.add_command("run_island", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_island运行完成")
        else:
            self.set_error("run_island运行出错!")

    def run_island_stat(self):
        out = self.work_dir + '/all.island.xls'
        cmd = '{} {}get_island.pl {} {} {}'.format(self.perl_path, self.perl_script,
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
        cmd = '{} {}get_phage_fasta.pl {} {} {} {} {}'.format(self.perl_path, self.perl_script, self.sample, "island",  self.scaf_seq,
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
        n = self.get_num(self.work_dir + '/all.island.xls')
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
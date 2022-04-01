#!/usr/bin/env python
# -*- coding: utf-8 -*-
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os, re
from mbio.packages.bacgenome.common import sum_stat
import unittest

class PanOrthomclAgent(Agent):
    """
    orthomcl做聚类生成同源基因
    author:qingchen.zhang
    """

    def __init__(self, parent):
        super(PanOrthomclAgent, self).__init__(parent)
        options = [
            {"name": "fasta_dir", "type": "infile", "format": "sequence.fasta_dir"},  #
            {"name": "pv_cutoff", "type": "float", "default": "1e-5"},  # P-Value or E-Value Cutoff
            {"name": "pi_cutoff", "type": "int", "default": 0, "max": 100, "min": 0},  # Percent Identity Cutoff
            {"name": "inflation", "type": "float", "default": 1.5},  # Markov Inflation Index
            {"name": "pmatch_cutoff", "type": "int", "default": 0, "max": 100, "min": 0},  # Percent Match Cutoff
            {"name": "result", "type": "outfile", "format": "sequence.profile_table"},
        ]
        self.add_option(options)
        self._memory_increase_step = 50  # 每次重运行增加内存50G by qingchen.zhang

    def check_options(self):
        """
        检测参数是否正确
        """
        if not self.option('fasta_dir').is_set:
            raise OptionError("设置fasta文件目录！")
        return True

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 4
        infile_path = self.option("infile_dir").prop['path']
        number_list = os.listdir(infile_path)
        number = len(number_list) / 2 ##几十个样本就用几十G
        total_memory = number * 1 + 10
        self._memory = '{}G'.format(total_memory)

    def end(self):
        super(PanOrthomclAgent, self).end()


class PanOrthomclTool(Tool):
    """
    version 1.0
    """
    def __init__(self, config):
        super(PanOrthomclTool, self).__init__(config)
        self.perl_path = "/program/perl-5.24.0/bin/perl"
        self.perl_script = self.config.PACKAGE_DIR + "/bacgenome/orthomcl.pl"
        self.orth_stat = self.config.PACKAGE_DIR + "/bac_comp_genome/orthomcl_stat.pl"
        self.blast = self.config.SOFTWARE_DIR + "/bioinfo/align/genBlast_v138_linux_x86_64/"
        self.mcl = self.config.SOFTWARE_DIR + "/bioinfo/Genomic/Sofware/mcl/bin"
        self.lib1 = self.config.SOFTWARE_DIR + "/program/perl-5.24.0/lib/5.24.0/BioPerl-1.6.924/"
        self.lib2 = self.config.SOFTWARE_DIR + "/bioinfo/Genomic/Script"
        self.set_environ(PATH=self.blast)
        self.set_environ(PATH=self.mcl)
        self.set_environ(PERL5LIB=self.lib1)
        self.set_environ(PERL5LIB=self.lib2)
        self.out = os.path.join(self.work_dir, 'orthomcl_out')

    def run_orthomcl(self):
        """
        orthomcl软件进行聚类
        """
        faa_dir = self.option("fasta_dir").prop["path"]
        cmd = '{} {} {} {} {} {} {} {}'.format(self.perl_path, self.perl_script, faa_dir,
                                                                            self.option('pv_cutoff'),
                                                                            self.option('pi_cutoff'),
                                                                            self.option('inflation'),
                                                                            self.option('pmatch_cutoff'), self.work_dir+"/")
        self.logger.info(cmd)
        command = self.add_command("orthomcl", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_orthomcl运行完成")
        else:
            self.set_error("run_orthomcl运行出错!")

    def creat_faa_list(self):
        """
        整理faa_list结果作为统计的输入
        """
        self.faa_list = os.path.join(self.work_dir,"faa_list")
        faa_dir = self.option("fasta_dir").prop["path"]
        all_files = os.listdir(faa_dir)
        with open(self.faa_list,"w") as f:
            for each in all_files:
                sample = each.split(".faa")[0]
                sample_path = os.path.join(faa_dir,each)
                f.write(sample + "\t" + sample_path + "\n")

    def run_stat(self):
        """
        开始对orthomcl的输出结果进行整理
        """
        self.logger.info("正在对结果进行统计和生成代表序列")
        self.mcl_file = os.path.join(self.work_dir, "all_orthomcl.out")
        if not os.path.exists(self.mcl_file):
            self.set_error("all_orthomcl.out creat failed")
        cmd = '{} {} {} {} {}'.format(self.perl_path,  self.orth_stat, self.mcl_file, self.faa_list, self.out)
        self.logger.info(cmd)
        command = self.add_command("run_stat", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_stat运行完成")
        else:
            self.set_error("run_stat运行出错!")

    def set_output(self):
        """
        设置结果文件目录
        """
        work_result = os.path.join(self.out, 'homologues_cluster.xls')
        result = os.path.join(self.output_dir, 'homologues_cluster.xls')
        if os.path.exists(result):
            os.remove(result)
        os.link(work_result, result)
        self.end()

    def run(self):
        """
        运行
        """
        super(PanOrthomclTool, self).run()
        self.creat_faa_list()
        self.run_orthomcl()
        self.run_stat()
        self.set_output()

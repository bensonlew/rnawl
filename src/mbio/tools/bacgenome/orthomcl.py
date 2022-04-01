#!/usr/bin/env python
# -*- coding: utf-8 -*-
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os, re
from mbio.packages.bacgenome.common import sum_stat
import unittest

class OrthomclAgent(Agent):
    """
    同源基因生成
    version 1.0
    author: ysh
    last_modify: 2019.04.16
    """

    def __init__(self, parent):
        super(OrthomclAgent, self).__init__(parent)
        options = [
            {"name": "fasta_dir", "type": "infile", "format": "sequence.fasta_dir"},  #
            {"name": "pv_cutoff", "type": "float", "default": "1e-5"},  # P-Value or E-Value Cutoff
            {"name": "pi_cutoff", "type": "int", "default": 0, "max": 100, "min": 0},  # Percent Identity Cutoff
            {"name": "inflation", "type": "float", "default": 1.5},  # Markov Inflation Index
            {"name": "pmatch_cutoff", "type": "int", "default": 0, "max": 100, "min": 0},  # Percent Match Cutoff
            {"name": "result", "type": "outfile", "format": "sequence.profile_table"},
        ]
        self.add_option(options)

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
        faa_dir = self.option("fasta_dir").prop["path"]
        all_files = os.listdir(faa_dir)
        self._memory = str(5 + len(all_files)*2) + "G"

    def end(self):
        super(OrthomclAgent, self).end()


class OrthomclTool(Tool):
    """
    version 1.0
    """
    def __init__(self, config):
        super(OrthomclTool, self).__init__(config)
        self.perl_path = "/program/perl-5.24.0/bin/perl"
        self.perl_script = self.config.PACKAGE_DIR + "/bacgenome/orthomcl.pl"
        self.orth_stat = self.config.PACKAGE_DIR + "/bacgenome/orthomcl_stat.pl"
        self.blast = self.config.SOFTWARE_DIR + "/bioinfo/align/genBlast_v138_linux_x86_64/"
        self.mcl = self.config.SOFTWARE_DIR + "/bioinfo/Genomic/Sofware/mcl/bin"
        self.lib1 = self.config.SOFTWARE_DIR + "/program/perl-5.24.0/lib/5.24.0/BioPerl-1.6.924/"
        self.lib2 = self.config.SOFTWARE_DIR + "/bioinfo/Genomic/Script"
        self.set_environ(PATH=self.blast)
        self.set_environ(PATH=self.mcl)
        self.set_environ(PERL5LIB=self.lib1)
        self.set_environ(PERL5LIB=self.lib2)

    def run_orthomcl(self):
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
        self.faa_list = os.path.join(self.work_dir,"faa_list")
        faa_dir = self.option("fasta_dir").prop["path"]
        all_files = os.listdir(faa_dir)
        with open(self.faa_list,"w") as f:
            for each in all_files:
                sample = each.split(".faa")[0]
                sample_path = os.path.join(faa_dir,each)
                f.write(sample + "\t" + sample_path + "\n")


    def run_stat(self):
        self.mcl_file = os.path.join(self.work_dir, "all_orthomcl.out")
        if not os.path.exists(self.mcl_file):
            self.set_error("all_orthomcl.out creat failed")
        cmd = '{} {} {} {} {}'.format(self.perl_path,  self.orth_stat, self.mcl_file, self.faa_list, self.output_dir)
        self.logger.info(cmd)
        command = self.add_command("run_stat", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_stat运行完成")
        else:
            self.set_error("run_stat运行出错!")

    def set_output(self):
        self.end()

    def run(self):
        """
        运行
        """
        super(OrthomclTool, self).run()
        self.creat_faa_list()
        self.run_orthomcl()
        self.run_stat()
        self.set_output()

    def get_num(self, file):
        with open(file, 'r') as f:
            lines = f.readlines()
            num = len(lines)
        return num

class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            #"id": "CreatTable" + str(random.randint(1, 10000)),
            "id": "Orthomcl",
            "type": "tool",
            "name": "bacgenome.orthomcl",
            "instant": False,
            "options": dict(
                fasta_dir="/mnt/ilustre/users/sanger-dev/sg-users/yuanshaohua/measurement/orthomcl/faa_dir",
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
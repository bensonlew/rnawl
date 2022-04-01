# -*- coding: utf-8 -*-
# __author__ = 'gao.hao'
# version 1.0
# last_modify: 2021.01.05

import os
import re
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import subprocess


class GenomeHousekeepingAgent(Agent):
    """
    prodigal 进行基因预测,再进行看家基因的
    """

    def __init__(self, parent):
        super(GenomeHousekeepingAgent, self).__init__(parent)
        options = [
            {"name": "input_genome", "type": "infile", "format": "sequence.fasta"},  # 组装拼接好的scaffold文件
            {"name": "sample", "type": "string", "default":"out"},
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("input_genome").is_set:
            raise OptionError("必须设置参数input_genome")

    def set_resource(self):
        self._cpu = 2
        self._memory = '10G'

    def end(self):
        super(GenomeHousekeepingAgent, self).end()


class GenomeHousekeepingTool(Tool):
    def __init__(self, config):
        super(GenomeHousekeepingTool, self).__init__(config)
        self.genome_fasta = self.option("input_genome").prop['path']
        self.prodigal_path = self.config.SOFTWARE_DIR +"/bioinfo/metaGenomic/Prodigal-2.6.3/prodigal"
        self.perl_path = "/program/perl-5.24.0/bin/perl"
        self.diamond = "/bioinfo/align/diamond-0.9.11/diamond"
        self.perl_script = self.config.PACKAGE_DIR + "/toolapps/"
        self.core_gene = self.config.SOFTWARE_DIR + "/database/bacgenome/House-keeping_gene/Core_gene"
        self.all_core_gene = self.config.SOFTWARE_DIR + "/database/bacgenome/House-keeping_gene/all_Coregene"
        self.all_core_gene_fa = self.config.SOFTWARE_DIR + "/database/bacgenome/House-keeping_gene/Core_gene.order.faa"
        self.list = self.config.SOFTWARE_DIR + "/database/bacgenome/House-keeping_gene/cor_gene.list"
        self.index = self.work_dir + "/ref"

    def run_prodigal(self):
        cmd = "{0} -i {2}  -o {1}.gff  -f gff -a {1}.faa -d {1}.ffn ".format(self.prodigal_path,self.option('sample'),self.genome_fasta)
        try:
            self.logger.info(cmd)
            subprocess.check_output(cmd, shell=True)
            self.logger.info("prodigal运行完成")
        except subprocess.CalledProcessError:
            self.set_error("prodigal运行出错")

    def run_index(self):
        cmd = "{} makedb --in {} -d {}".format(self.diamond, self.work_dir+'/'+self.option('sample')+'.faa', self.index)
        self.logger.info(cmd)
        self.logger.info("开始运行run_index")
        command = self.add_command("run_index", cmd)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("运行run_index完成")
        else:
            self.set_error("运行run_index运行出错!")

    def run_get_gene(self):
        cmd = "{} blastp -d {} -q {} -o {} --max-target-seqs 1".format(self.diamond, self.index,self.all_core_gene_fa,self.work_dir +'/' + self.option('sample') + ".matches.m8")
        self.logger.info(cmd)
        self.logger.info("开始运行run_get_gene")
        command = self.add_command("run_get_gene", cmd)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("运行run_get_gene完成")
        else:
            self.set_error("运行run_get_gene运行出错!", code="31402001")

    def run_get_seq(self):
        cmd = "{} {}combin_cor_gene.pl {} {} {} {} {}".format(self.perl_path, self.perl_script,self.option('sample'),self.work_dir +'/' + self.option('sample') + ".matches.m8" , self.work_dir+'/'+self.option('sample')+'.faa', self.list,self.output_dir + '/' + self.option('sample') + '.core_gene.fa')
        self.logger.info(cmd)
        self.logger.info("开始运行run_get_seq")
        command = self.add_command("run_get_seq", cmd)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("运行run_get_seq完成")
        else:
            self.set_error("运行run_get_seq运行出错!")

    def set_output(self):
        if os.path.exists(self.output_dir +'/' + self.option('sample') + ".matches.m8"):
            os.remove(self.output_dir +'/' + self.option('sample') + ".matches.m8")
        os.link(self.work_dir +'/' + self.option('sample') + ".matches.m8", self.output_dir +'/' + self.option('sample') + ".matches.m8")

    def run(self):
        super(GenomeHousekeepingTool, self).run()
        self.run_prodigal()
        self.run_index()
        self.run_get_gene()
        self.run_get_seq()
        self.set_output()
        self.end()
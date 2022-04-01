#!/usr/bin/env python
# -*- coding: utf-8 -*-
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os,re
from mbio.packages.metagbin.common_function import link_dir


class GetCoreGeneAgent(Agent):
    """
    用于获取看家基因进化树序列
    version 1.0
    author: gaohao
    last_modify: 2018.04.21
    """

    def __init__(self, parent):
        super(GetCoreGeneAgent, self).__init__(parent)
        options = [
            {"name": "seq_faa", "type": "infile", "format": "sequence.fasta"},
            {"name": "seq_gff", "type": "infile", "format": "gene_structure.gff3,metagbin.file_table"},
            {"name": "sample_name", "type": "string"},
            {"name": "method", "type": "string", "default": "bins"} ##bins or genomes
        ]
        self.add_option(options)
        self.list =[]

    def check_options(self):
        """
        检测参数是否正确
        """
        if not self.option('seq_faa').is_set:
            raise OptionError("请设置基因组基因蛋白序列不存在！", code="31402001")

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 5
        self._memory = '40G'

    def end(self):
        super(GetCoreGeneAgent, self).end()

class GetCoreGeneTool(Tool):
    """
    version 1.0
    """
    def __init__(self, config):
        super(GetCoreGeneTool, self).__init__(config)
        self.gff = self.option('seq_gff').prop['path']
        self.gene =self.option('seq_faa').prop['path']
        self.perl_path = "/program/perl-5.24.0/bin/perl"
        self.perl_script = self.config.PACKAGE_DIR + "/bac_comp_genome/"
        self.diamond = "/bioinfo/align/diamond-0.8.35/diamond"
        self.core_gene = self.config.SOFTWARE_DIR + "/database/bacgenome/House-keeping_gene/Core_gene"
        self.all_core_gene = self.config.SOFTWARE_DIR + "/database/bacgenome/House-keeping_gene/all_Coregene"
        self.all_core_gene_fa = self.config.SOFTWARE_DIR + "/database/bacgenome/House-keeping_gene/Core_gene.order.faa"
        self.list = self.config.SOFTWARE_DIR + "/database/bacgenome/House-keeping_gene/cor_gene.list"
        self.index = self.work_dir + "/ref"

    def run_index(self):
        cmd = "{} makedb --in {} -d {}".format(self.diamond, self.gene, self.index)
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
        if os.path.exists(self.work_dir + '/all.fasta'):
            os.remove(self.work_dir + '/all.fasta')
        cmd = "{} blastp -d {} -q {} -o {} --max-target-seqs 1".format(self.diamond, self.index,self.all_core_gene_fa,self.work_dir +'/' + self.option('sample_name') + ".matches.m8")
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
        cmd = "{} {}combin_cor_gene.pl {} {} {} {} {}".format(self.perl_path, self.perl_script,self.option('sample_name'),self.work_dir +'/' + self.option('sample_name') + ".matches.m8" , self.gene, self.list,self.work_dir + '/' + self.option('sample_name') + '.cor_gene.fa')
        self.logger.info(cmd)
        self.logger.info("开始运行run_get_seq")
        command = self.add_command("run_get_seq", cmd)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("运行run_get_seq完成")
        else:
            self.set_error("运行run_get_seq运行出错!")

    def run_getresult(self):
        if self.option("method") in ["bins"]:
            cmd = "{} {}tiqu_bins.pl {} {} {} {} {}".format(self.perl_path, self.perl_script,
                                                              self.option('sample_name'),
                                                              self.work_dir + '/' + self.option(
                                                                  'sample_name') + ".matches.m8", self.gene, self.gff,
                                                              self.work_dir + '/' + self.option(
                                                                  'sample_name') + '.result.xls')
        else:
            cmd = "{} {}tiqu_genomes.pl {} {} {} {} {}".format(self.perl_path, self.perl_script,
                                                            self.option('sample_name'),
                                                            self.work_dir + '/' + self.option(
                                                                'sample_name') + ".matches.m8", self.gene, self.gff,
                                                            self.work_dir + '/' + self.option(
                                                                'sample_name') + '.result.xls')
        self.logger.info(cmd)
        self.logger.info("开始运行run_getresult")
        command = self.add_command("run_getresult", cmd)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("运行run_getresult完成")
        else:
            self.set_error("运行run_getresult运行出错!")

    def set_output(self):
        for i in [self.option('sample_name') + ".result.xls",self.option('sample_name') + ".cor_gene.fa"]:
            if os.path.exists(self.output_dir + "/" + i):
                os.remove(self.output_dir + "/" + i)
            os.link(self.work_dir + "/" + i, self.output_dir + "/" + i)

    def run(self):
        """
        运行
        """
        super(GetCoreGeneTool, self).run()
        self.run_index()
        self.run_get_gene()
        self.run_get_seq()
        self.run_getresult()
        self.set_output()
        self.end()
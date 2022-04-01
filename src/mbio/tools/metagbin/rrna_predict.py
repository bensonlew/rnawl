# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
# version 1.0
# last_modify: 20119.01.14

import os
import re
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import subprocess


class RrnaPredictAgent(Agent):
    """
    Barrnap 进行宏基因组binning结果rRNA预测
    """

    def __init__(self, parent):
        super(RrnaPredictAgent, self).__init__(parent)
        options = [
            {"name": "bin_dir", "type": "infile", "format": "sequence.fasta_dir"},  # bin的文件目录
            {"name": "input_genome", "type": "infile", "format": "sequence.fasta"},  # 组装拼接好的scaffold文件
            {"name": "kingdom", "type": "string", "default": "bac"},  # genome type Kingdom: arc mito bac euk
            {"name": "lencutoff", "type": "string", "default": "0.8"},
            # Proportional length threshold to label as partial，小于该阈值预测结果会加上(partial)
            {"name": "reject", "type": "string", "default": "0.5"},
            {"name": "metabat_depth", "type": "outfile", "format": "sequence.profile_table"},
            # Proportional length threshold to reject prediction， 只有当大于该阈值的结果才会保留
            {"name": "evalue", "type": "float", "default": 1e-06},  # Similarity e-value cut-off
            {"name": "seq", "type": "outfile", "format": "sequence.fasta"},  # 预测rRNA序列
            {"name": "rna_gff", "type": "outfile", "format": "gene_structure.gff3"},  # 预测生成的GFF格式,并增加两列位置信息，用与绘制圈图
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("input_genome").is_set:
            raise OptionError("必须设置参数input_genome")

    def set_resource(self):
        self._cpu = 4
        self._memory = '20G'

    def end(self):
        super(RrnaPredictAgent, self).end()


class RrnaPredictTool(Tool):
    def __init__(self, config):
        super(RrnaPredictTool, self).__init__(config)
        self.genome_fasta = self.option("input_genome").prop['path']
        self.perl_path = "/program/perl-5.24.0/bin/perl"
        self.perl_script = self.config.PACKAGE_DIR + "/sequence/scripts/"
        self.perl_script2 = self.config.PACKAGE_DIR + "/metagbin/"
        self.barrnap_path = self.config.SOFTWARE_DIR +"/bioinfo/Genomic/Sofware/barrnap-0.8/bin/"
        self.scf_sum = self.work_dir + '/all_bins.scf.xls'
        self.s16_sum = self.work_dir + '/all_bins.16s.xls'

    def run_barrnap(self):
        cmd = "{}barrnap --threads 2 --kingdom {} --quiet {} --lencutoff {} --reject {} --evalue {} > {}".format(
            self.barrnap_path,
            self.option("kingdom"), self.genome_fasta,
            self.option("lencutoff"), self.option("reject"),
            self.option("evalue"), "bins.gff")
        try:
            self.logger.info(cmd)
            subprocess.check_output(cmd, shell=True)
            self.logger.info("barrnap运行完成")
        except subprocess.CalledProcessError:
            self.set_error("barrnap运行出错")

    def run_rrna_extract(self):
        gff_file = self.work_dir + "/bins.gff"
        cmd = '{} {}noncrna_fnn_extract.pl {} {} {} {}'.format(self.perl_path, self.perl_script, self.genome_fasta,
                                                               gff_file, "rrna", self.work_dir)
        command = self.add_command("noncrna_fnn_extract", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("提取fnn和整理GFF文件运行完成")
        else:
            self.set_error("提取fnn和整理GFF文件运行出错!")

    def run_summary(self):
        bin_dir =self.option("bin_dir").prop['path']
        metabat_depth = self.option("metabat_depth").prop['path']
        cmd = '{} {}bin_scaffold_sum.pl {} {} {} {}'.format(self.perl_path, self.perl_script2,bin_dir, self.genome_fasta,metabat_depth,self.scf_sum)
        command = self.add_command("run_summary", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("提取run_summary信息运行完成")
        else:
            self.set_error("提取frun_summary信息运行出错!")

    def run_rrna_sum(self):
        bin_dir =self.option("bin_dir").prop['path']
        cmd = '{} {}scaffold_16s.pl {} {} {} {}'.format(self.perl_path, self.perl_script2,bin_dir, self.genome_fasta,self.work_dir + "/" + "bins.gff",self.s16_sum)
        command = self.add_command("run_rrna_sum", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("提取run_rrna_sum信息运行完成")
        else:
            self.set_error("提取run_rrna_sum信息运行出错!")

    def set_output(self):
        for i in ['bins.rRNA.gff','bins.rRNA.fnn','all_bins.scf.xls','all_bins.16s.xls']:
            if os.path.exists(self.output_dir + '/' + i):
                os.remove(self.output_dir + '/' + i)
            os.link(self.work_dir + '/' + i,self.output_dir + '/' +  i)

    def run(self):
        super(RrnaPredictTool, self).run()
        self.run_barrnap()
        self.run_rrna_extract()
        self.run_summary()
        self.run_rrna_sum()
        self.set_output()
        self.end()

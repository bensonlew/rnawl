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

class BinRrnaAgent(Agent):
    """
    Barrnap 进行bin结果rRNA预测
    """

    def __init__(self, parent):
        super(BinRrnaAgent, self).__init__(parent)
        options = [
            {"name": "bin_fa", "type": "infile", "format": "sequence.fasta"},  # bin的fasta文件
            {"name": "kingdom", "type": "string", "default": "bac"},  # genome type Kingdom: arc mito bac euk
            {"name": "lencutoff", "type": "string", "default": "0.8"},
            # Proportional length threshold to label as partial，小于该阈值预测结果会加上(partial)
            {"name": "reject", "type": "string", "default": "0.5"},
            {"name": "evalue", "type": "float", "default": 1e-06},  # Similarity e-value cut-off
            {"name": "seq", "type": "outfile", "format": "sequence.fasta"},  # 预测rRNA序列
            {"name": "rna_gff", "type": "outfile", "format": "gene_structure.gff3"},  # 预测生成的GFF格式,并增加两列位置信息，用与绘制圈图
        ]
        self.add_option(options)
        self.type = ''

    def check_options(self):
        if not self.option("bin_fa").is_set:
            raise OptionError("必须设置参数bin_fa")
        if not self.option("kingdom"):
            raise OptionError("必须设置参数kingdom")

    def set_resource(self):
        self._cpu = 2
        self._memory = '2G'

    def end(self):
        super(BinRrnaAgent, self).end()


class BinRrnaTool(Tool):
    def __init__(self, config):
        super(BinRrnaTool, self).__init__(config)
        self.genome_fasta = self.option("bin_fa").prop['path']
        self.sample_name = os.path.basename(self.option("bin_fa").prop['path']).split(".")[0]
        self.perl_path = "/program/perl-5.24.0/bin/perl"
        self.perl_script = self.config.PACKAGE_DIR + "/sequence/scripts/"
        self.perl_script2 = self.config.PACKAGE_DIR + "/metagbin/"
        self.barrnap_path = self.config.SOFTWARE_DIR +"/bioinfo/Genomic/Sofware/barrnap-0.8/bin/"
        self.scf_sum = self.work_dir + '/all_bins.scf.xls'
        self.s16_sum = self.work_dir + '/all_bins.16s.xls'
        if self.option("kingdom") in ['Bacteria']:
            self.kingdom = 'bac'
        elif self.option("kingdom") in ['Archaea']:
            self.kingdom = "arc"

    def run_barrnap(self):
        if os.path.exists(self.work_dir + '/' + self.sample_name + ".gff"):
            os.remove(self.work_dir + '/' + self.sample_name + ".gff")
        cmd = "{}barrnap --threads 2 --kingdom {} --quiet {} --lencutoff {} --reject {} --evalue {} > {}".format(
            self.barrnap_path,self.kingdom, self.genome_fasta,
            self.option("lencutoff"), self.option("reject"),
            self.option("evalue"), self.sample_name + ".gff")
        try:
            self.logger.info(cmd)
            subprocess.check_output(cmd, shell=True)
            self.logger.info("barrnap运行完成")
        except subprocess.CalledProcessError:
            self.set_error("barrnap运行出错")

    def run_rrna_extract(self):
        gff_file = self.work_dir + "/" +self.sample_name + ".gff"
        cmd = '{} {}noncrna_fnn_extract.pl {} {} {} {}'.format(self.perl_path, self.perl_script, self.genome_fasta,
                                                               gff_file, "rrna", self.work_dir)
        command = self.add_command("noncrna_fnn_extract", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.type = self.get_s16(self.work_dir + "/" + self.sample_name + ".rRNA.gff")
            self.logger.info("提取fnn和整理GFF文件运行完成")
        else:
            self.set_error("提取fnn和整理GFF文件运行出错!")

    def run_rrna_sum(self):
        bin_dir =os.path.dirname(self.option("bin_fa").prop['path'])
        cmd = '{} {}scaffold_16s.pl {} {} {} {}'.format(self.perl_path, self.perl_script2,bin_dir, self.genome_fasta,self.work_dir + "/" + self.sample_name + ".gff",self.s16_sum)
        command = self.add_command("run_rrna_sum", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("提取run_rrna_sum信息运行完成")
        else:
            self.set_error("提取run_rrna_sum信息运行出错!")

    def set_output(self):
        for i in [self.sample_name + '.rRNA.gff',self.sample_name + '.rRNA.fna','all_bins.16s.xls']:
            if os.path.exists(self.output_dir + '/' + i):
                os.remove(self.output_dir + '/' + i)
            if os.path.exists(self.work_dir + '/' + i):
                os.link(self.work_dir + '/' + i,self.output_dir + '/' +  i)

    def get_size(self,file):
        with open(file, "r") as f:
            lines =f.readlines()
            size =len(lines)
        return size

    def get_s16(self,file):
        type = 'false'
        with open (file,'r') as f:
            lines = f.readlines()
            for line in lines[1:]:
                if re.search(r'ame=16S_rRNA',line):
                    type = "true"
            return type

    def run(self):
        super(BinRrnaTool, self).run()
        self.run_barrnap()
        if self.get_size(self.work_dir + "/" +self.sample_name + ".gff") >=2:
            self.run_rrna_extract()
            if self.type in ['true']:
                self.run_rrna_sum()
                self.set_output()
                self.end()
            else:
                self.end()
        else:
            self.end()


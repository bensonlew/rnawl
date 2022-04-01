# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
# version 1.0
# last_modify: 2021.01.05

import os
import re
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import subprocess


class GenomeRrnaAgent(Agent):
    """
    Barrnap 进行rRNA预测
    """

    def __init__(self, parent):
        super(GenomeRrnaAgent, self).__init__(parent)
        options = [
            {"name": "input_genome", "type": "infile", "format": "sequence.fasta"},  # 组装拼接好的scaffold文件
            {"name": "kingdom", "type": "string", "default": "bac"},  # genome type Kingdom: arc mito bac euk
            {"name": "lencutoff", "type": "string", "default": "0.8"},
            # Proportional length threshold to label as partial，小于该阈值预测结果会加上(partial)
            {"name": "reject", "type": "string", "default": "0.5"},
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
        self._cpu = 2
        self._memory = '2G'

    def end(self):
        super(GenomeRrnaAgent, self).end()


class GenomeRrnaTool(Tool):
    def __init__(self, config):
        super(GenomeRrnaTool, self).__init__(config)
        self.genome_fasta = self.option("input_genome").prop['path']
        self.sample_name = (self.genome_fasta).split('/')[-1].rsplit('.')[0]
        self.perl_path = "/program/perl-5.24.0/bin/perl"
        self.perl_script = self.config.PACKAGE_DIR + "/toolapps/"
        self.barrnap_path = self.config.SOFTWARE_DIR +"/bioinfo/Genomic/Sofware/barrnap-0.8/bin/"

    def run_barrnap(self):
        cmd = "{}barrnap --threads 2 --kingdom {} --quiet {} --lencutoff {} --reject {} --evalue {} > {}".format(
            self.barrnap_path,
            self.option("kingdom"), self.genome_fasta,
            self.option("lencutoff"), self.option("reject"),
            self.option("evalue"), self.sample_name + ".gff")
        try:
            self.logger.info(cmd)
            subprocess.check_output(cmd, shell=True)
            self.logger.info("barrnap运行完成")
            self.run_rrna_extract()
        except subprocess.CalledProcessError:
            self.set_error("barrnap运行出错", code="33300201")

    def run_rrna_extract(self):
        gff_file = self.work_dir + "/" +self.sample_name + ".gff"
        cmd = '{} {}noncrna_fnn_extract.pl {} {} {} {} {} {}'.format(self.perl_path, self.perl_script, self.genome_fasta,
                                                               gff_file, "rrna", self.output_dir, self.sample_name, "gene")
        command = self.add_command("noncrna_fnn_extract", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("提取fnn和整理GFF文件运行完成")
            self.set_output()
        else:
            self.set_error("提取fnn和整理GFF文件运行出错!", code="33300202")

    def set_output(self):
        nul_path = self.output_dir + "/" + self.sample_name + "_16S.ffn"
        if os.path.getsize(nul_path):
            self.option('seq', nul_path)

    def run(self):
        super(GenomeRrnaTool, self).run()
        self.run_barrnap()
        self.end()
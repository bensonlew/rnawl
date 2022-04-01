#-*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'@20190912

import os
import re
import shutil
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import subprocess


class BarrnapAgent(Agent):
    """
    细菌比较基因组 Barrnap 进行rRNA预测可以用gff进行提取序列
    1. 先用软件Barrnap进行预测；
    2. 根据预测结果提取序列和矫正
    """

    def __init__(self, parent):
        super(BarrnapAgent, self).__init__(parent)
        options = [
            {"name": "input_genome", "type": "infile", "format": "sequence.fasta"},  # 组装拼接好的scaffold文件
            {"name": "kingdom", "type": "string", "default": "bac"},  # genome type Kingdom: arc mito bac euk
            {"name": "lencutoff", "type": "string", "default": "0.8"},
            # Proportional length threshold to label as partial，小于该阈值预测结果会加上(partial)
            {"name": "reject", "type": "string", "default": "0.5"},
            # Proportional length threshold to reject prediction， 只有当大于该阈值的结果才会保留
            {"name": "evalue", "type": "float", "default": 1e-06},  # Similarity e-value cut-off
            {"name": "genome_name", "type": "string"},  # 样品下的基因组名称NZ_CM000741.1或者GCF_000003955.1
            {"name": "genome_type", "type": "string"},  # 基因组属于什么类型，complete，draft,chromosome
            {"name": "gene_tag", "type": "string"}, # 基因前缀为多少，主要是用于提取序列用和整理成统一的序列用的
            {"name": "seq", "type": "outfile", "format": "sequence.fasta"},  # 预测rRNA序列
            {"name": "gff", "type": "outfile", "format": "gene_structure.gff3"},  # 预测生成的GFF格式,并增加两列位置信息，用与绘制圈图
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("genome_name"):
            raise OptionError("必须设置参数genome_name")
        if not self.option("input_genome").is_set:
            raise OptionError("必须设置参数input_genome")
        if not self.option("gene_tag"):
            raise OptionError("必须设置参数gene_tag")

    def set_resource(self):
        self._cpu = 2
        self._memory = '2G'

    def end(self):
        super(BarrnapAgent, self).end()


class BarrnapTool(Tool):
    def __init__(self, config):
        super(BarrnapTool, self).__init__(config)
        self.genome_fasta = self.option("input_genome").prop['path']
        self.sample_name = self.option('genome_name')
        self.perl_path = "/program/perl-5.24.0/bin/perl"
        self.perl_script = self.config.PACKAGE_DIR + "/sequence/scripts/"
        self.barrnap_path = self.config.SOFTWARE_DIR +"/bioinfo/Genomic/Sofware/barrnap-0.8/bin/"
        self.package_path = self.config.PACKAGE_DIR + '/bac_comp_genome/'

    def run_barrnap(self):
        """
        根据序列用barrnap进行rRNA预测
        :return:
        """
        cmd = "{}barrnap --threads 2 --kingdom {} --quiet {} --lencutoff {} --reject {} --evalue {} > {}".format(
            self.barrnap_path,
            self.option("kingdom"), self.genome_fasta,
            self.option("lencutoff"), self.option("reject"),
            self.option("evalue"), self.sample_name + ".gff")
        try:
            self.logger.info(cmd)
            subprocess.check_output(cmd, shell=True)
            self.logger.info("barrnap运行完成")
        except subprocess.CalledProcessError:
            self.set_error("barrnap运行出错")

    def run_rrna_extract(self):
        """
        给原始的gff文件序列名称和location
        :return:
        """
        gff_file = self.work_dir + "/" +self.sample_name + ".gff"
        cmd = '{} {}noncrna_fnn_extract.pl {} {} {} {} {} {}'.format(self.perl_path, self.package_path, self.genome_fasta,
                                                               gff_file, "rrna", self.output_dir, self.sample_name, self.option('gene_tag'))
        self.logger.info(cmd)
        command = self.add_command("noncrna_fnn_extract", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("提取fnn和整理GFF文件运行完成")
            self.set_output()
        else:
            self.set_error("提取fnn和整理GFF文件运行出错!")

    def set_output(self):
        """
        设置结果文件目录
        :return:
        """
        self.logger.info("正在设置结果文件目录")
        nul_path = self.output_dir + "/" + self.sample_name + "_rRNA.fna"
        gff = self.output_dir + "/" + self.sample_name + "_rRNA.gff"
        if os.path.getsize(nul_path):
            self.option('seq', nul_path)
        self.option('gff', gff)
        self.logger.info('设置结果文件目录完成')

    def run(self):
        self.logger.info('开始运行啦')
        super(BarrnapTool, self).run()
        self.run_barrnap()
        gff_file = self.work_dir + "/" +self.sample_name + ".gff"
        with open(gff_file, 'r') as f:
            lines = f.readlines()
            line_num = len(lines)
            if line_num >= 2:
                self.run_rrna_extract()
                self.end()
            else:
                self.end()

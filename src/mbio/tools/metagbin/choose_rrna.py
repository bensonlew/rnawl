# -*- coding: utf-8 -*-
# __author__ = "qingchen.zhang"@20190127

from biocluster.agent import Agent
from biocluster.tool import Tool
import os
import subprocess
from biocluster.core.exceptions import OptionError
import shutil
import re


class ChooseRrnaAgent(Agent):
    """
    从rRNA预测的序列中挑选出16sRNA
    """
    def __init__(self, parent):
        super(ChooseRrnaAgent, self).__init__(parent)
        options = [
            {"name": "input_genome", "type": "infile", "format": "sequence.fasta"},  # 预测rRNA序列
            {"name": "rrna_gff", "type": "infile", "format": "gene_structure.gff3"},  # 预测生成的GFF格式,并增加两列位置信息，用与绘制圈图
            {"name": "sample", "type": "string"},   #基因组名称
            {"name": "seq", "type": "outfile", "format": "sequence.fasta"},  # 预测筛选出来的16S_rRNA序列
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("input_genome").is_set:
            raise OptionError("必须设置参数input_genome")
        if not self.option("rrna_gff").is_set:
            raise OptionError("必须设置参数rrna_gff")

    def set_resource(self):
        self._cpu = 2
        self._memory = '2G'

    def end(self):
        super(ChooseRrnaAgent, self).end()


class ChooseRrnaTool(Tool):
    def __init__(self, config):
        super(ChooseRrnaTool, self).__init__(config)
        self.perl_path = self.config.SOFTWARE_DIR + '/program/perl-5.24.0/bin/perl'
        self.script = self.config.PACKAGE_DIR + '/metagbin/'

    def run_creat_list(self):
        """
        从gff文件中挑选出16s_RNA
        :return:
        """
        stat_path = self.option('rrna_gff').prop['path']
        list_path = self.work_dir + '/16S_rRNA.txt'
        with open(stat_path, 'r') as infile, open(list_path, 'w') as outfile:
            lines = infile.readlines()
            for line in lines[1:]:
                line = line.strip().split('\t')
                attributes = line[7]

                if re.search(r'16S_rRNA', attributes):
                    gene_id = line[0]
                    outfile.write('{} '.format(gene_id))

    def run_choose(self):
        """
        根据list文件，从fasta中挑出16s
        """
        rna_path = self.option('input_genome').prop['path']
        if os.path.exists(self.work_dir + '/16S_rRNA.txt') and os.path.getsize(self.work_dir + '/16S_rRNA.txt') != 0:
            list_path = self.work_dir + '/16S_rRNA.txt'
            cmd = '{} {}choose_seqs.pl -f {} -l {} -o {}'.format(self.perl_path, self.script, rna_path, list_path, self.output_dir +'/'+self.option('sample') + '-16S_rRNA.fa')
            try:
                self.logger.info(cmd)
                subprocess.check_output(cmd, shell=True)
                self.logger.info("16s提取成功")
            except subprocess.CalledProcessError:
                self.logger.error("未能成功预测出16s_rRNA")
        else:
            self.logger.error("未能成功的预测出16s_rRNA")

    def set_output(self):
        """
        设置结果目录
        """
        self.logger.info("正在生成结果目录")
        if os.path.exists(self.output_dir +'/'+self.option('sample') + '-16S_rRNA.fa'):
            self.logger.info("生成结果目录成功")
            self.option('seq', self.output_dir +'/'+self.option('sample') + '-16S_rRNA.fa')
        else:
            self.logger.info('未能成功预测出16s_rRNA')

    def run(self):
        super(ChooseRrnaTool, self).run()
        self.run_creat_list()
        self.run_choose()
        self.set_output()
        self.end()

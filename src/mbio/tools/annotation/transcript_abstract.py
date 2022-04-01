# -*- coding: utf-8 -*-
# __author__ = 'zengjing'

import os
import shutil
import re
import subprocess
from biocluster.core.exceptions import OptionError
from biocluster.agent import Agent
from biocluster.tool import Tool

class TranscriptAbstractAgent(Agent):
    """
    提取参考基因组最长序列，作为基因注释的输入文件
    author: zengjing
    last_modify: 2016.09.22
    """
    def __init__(self, parent):
        super(TranscriptAbstractAgent, self).__init__(parent)
        options = [
            {"name": "ref_genome_custom", "type": "infile", "format": "sequence.fasta"},  # 参考基因组fasta文件
            {"name": "ref_genome_gtf", "type": "infile", "format": "gene_structure.gtf"},  # 参考基因组gtf文件
            {"name": "ref_genome_gff", "type": "infile", "format": "gene_structure.gff3"},  # 参考基因组gff文件
            {"name": "query", "type": "outfile", "format": "sequence.fasta"},  # 输出做注释的转录本序列
            {"name": "gene_file", "type": "outfile", "format": "rna.gene_list"},  # 输出最长转录本
            {"name": "length_file", "type": "outfile", "format": "annotation.cog.cog_list"}  # 输出注释转录本序列的长度
        ]
        self.add_option(options)
        self.step.add_steps("Transcript")
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.Transcript.start()
        self.step.update()

    def step_end(self):
        self.step.Transcript.finish()
        self.step.update()

    def check_option(self):
        if not self.option("ref_genome_custom").is_set:
            raise OptionError("请设置参考基因组custom文件")
        if self.option("ref_genome_gtf").is_set or self.option("ref_genome_gff").is_set:
            pass
        else:
            raise OptionError("请设置参考基因组gtf文件或gff文件")

    def set_resource(self):
        self._cpu = 1
        self._memory = '5G'

    def end(self):
        super(TranscriptAbstractAgent, self).end()

class TranscriptAbstractTool(Tool):
    def __init__(self, config):
        super(TranscriptAbstractTool, self).__init__(config)
        self.gffread_path = "bioinfo/rna/cufflinks-2.2.1/"
        self.long_path = self.config.SOFTWARE_DIR + "/bioinfo/rna/scripts/"
        self.python_path = "program/Python/bin/"
        self.length = self.config.SOFTWARE_DIR + "/bioinfo/annotation/scripts/fastalength"

    def run_gffread(self):
        old_fasta = self.option("ref_genome_custom").prop["path"]  # 统一按自定义方式传参考基因组
        fasta = self.work_dir + "/" + os.path.basename(self.option("ref_genome_custom").prop["path"])
        if os.path.exists(fasta):
            os.remove(fasta)
        os.link(old_fasta, fasta)
        gtf = self.option("ref_genome_gtf").prop["path"]
        cmd = "{}gffread {} -g {} -w exons.fa".format(self.gffread_path, gtf, fasta)
        self.logger.info("开始运行cufflinks的gffread，合成、提取exons")
        command = self.add_command("gffread", cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("exons提取完成")
        else:
            self.set_error("运exons提取出错")

    def run_exons_length(self):
        """exons的长度"""
        exon_path = os.path.join(self.work_dir, "exons.fa")
        length_path = self.output_dir + "/exons_length.txt"
        cmd = "{} {} > {}".format(self.length, exon_path, length_path)
        self.logger.info("开始提取exons的长度")
        try:
            subprocess.check_output(cmd, shell=True)
            self.logger.info("提取exons的长度完成")
        except subprocess.CalledProcessError:
            self.set_error("提取exons的长度失败")
        self.option("length_file", length_path)

    def run_long_transcript(self):
        exon_path = os.path.join(self.work_dir, "exons.fa")
        cmd = "{}python {}annotation_longest.py -i {}".format(self.python_path, self.long_path, exon_path)
        self.logger.info("提取最长序列")
        command = self.add_command("the_longest", cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("运提取最长序列完成")
        else:
            self.set_error("提取最长序列出错")
        output1 = os.path.join(self.work_dir, "exons.fa")
        if os.path.exists(self.output_dir + "/exons.fa"):
            os.remove(self.output_dir + "/exons.fa")
        shutil.move(output1, self.output_dir + "/exons.fa")
        output2 = os.path.join(self.work_dir, "the_longest_exons.fa")
        if os.path.exists(self.output_dir + "/the_longest_exons.fa"):
            os.remove(self.output_dir + "/the_longest_exons.fa")
        shutil.move(output2, self.output_dir + "/the_longest_exons.fa")
        self.option('query', self.work_dir + '/output/exons.fa')

    def get_gene_list(self):
        output_path = self.work_dir + "/output/the_longest_exons.fa"
        gene_list_path = self.work_dir + '/output/gene_list.txt'
        gene_lists = []
        with open(output_path, 'rb') as f, open(gene_list_path, 'wb') as w:
            lines = f.readlines()
            for line in lines:
                m = re.match(r">(.+) gene=(.+)", line)
                if m:
                    trans_name = m.group(1)
                    if trans_name not in gene_lists:
                        w.write(trans_name + '\n')
                        gene_lists.append(trans_name)
                else:
                    n = re.match(r">(.+) transcript:(.+)", line)
                    if n:
                        trans_name = n.group(1)
                        if trans_name not in gene_lists:
                            w.write(trans_name + '\n')
                            gene_lists.append(trans_name)
        self.option('gene_file', gene_list_path)

    def run(self):
        super(TranscriptAbstractTool, self).run()
        self.run_gffread()
        self.run_exons_length()
        self.run_long_transcript()
        self.get_gene_list()
        self.end()

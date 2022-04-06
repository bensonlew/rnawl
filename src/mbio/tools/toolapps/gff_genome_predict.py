# -*- coding: utf-8 -*-
# __author__ = 'gao.hao'
# version 1.0
# last_modify: 2021.01.06

import os
import re
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import subprocess
from Bio import SeqIO


class GffGenomePredictAgent(Agent):
    """
    prodigal 进行基因预测
    """

    def __init__(self, parent):
        super(GffGenomePredictAgent, self).__init__(parent)
        options = [
            {"name": "input_genome", "type": "infile", "format": "sequence.fasta"},  # 组装拼接好的scaffold文件
            {"name": "gff", "type": "infile", "format": "gene_structure.gff3"},
            {"name": "sample", "type": "string","default":"out"},
            {"name": "s16", "type": "outfile", "format": "sequence.fasta"},
            {"name": "house", "type": "outfile", "format": "sequence.fasta"},
            {'name': 'stat', 'type': 'outfile', "format": "sequence.profile_table"},
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("input_genome").is_set:
            raise OptionError("必须设置参数input_genome")

    def set_resource(self):
        self._cpu = 2
        self._memory = '15G'

    def end(self):
        super(GffGenomePredictAgent, self).end()


class GffGenomePredictTool(Tool):
    def __init__(self, config):
        super(GffGenomePredictTool, self).__init__(config)
        self.genome_fasta = self.option("input_genome").prop['path']
        self.gff = self.option("gff").prop['path']
        self.prodigal_path = self.config.SOFTWARE_DIR +"/bioinfo/metaGenomic/Prodigal-2.6.3/prodigal"
        self.perl_path = "/program/perl-5.24.0/bin/perl"
        self.python_path = "/miniconda2/bin/python"
        self.perl_script = self.config.PACKAGE_DIR + "/toolapps/"
        self.diamond = "/bioinfo/align/diamond-0.8.35/diamond"
        self.transeq = "/bioinfo/seq/EMBOSS-6.6.0/emboss/transeq"
        self.LD_LIBRARY_PATH = self.config.SOFTWARE_DIR + "/bioinfo/seq/EMBOSS-6.6.0/lib"
        self.set_environ(LD_LIBRARY_PATH=self.LD_LIBRARY_PATH)
        self.core_gene = self.config.SOFTWARE_DIR + "/database/bacgenome/House-keeping_gene/Core_gene"
        self.all_core_gene = self.config.SOFTWARE_DIR + "/database/bacgenome/House-keeping_gene/all_Coregene"
        self.all_core_gene_fa = self.config.SOFTWARE_DIR + "/database/bacgenome/House-keeping_gene/Core_gene.order.faa"
        self.list = self.config.SOFTWARE_DIR + "/database/bacgenome/House-keeping_gene/cor_gene.list"
        self.index = self.work_dir + "/ref"

    def run_getseq(self):
        cmd = "{} {}get_seq_from_genome.py -s {}  -g {}  -p {}".format(self.python_path, self.perl_script, self.genome_fasta, self.gff, self.option('sample'))
        self.logger.info(cmd)
        self.logger.info("开始运行run_getseq")
        command = self.add_command("run_getseq", cmd)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("运行run_getseq完成")
        else:
            self.set_error("运行run_getseq运行出错!")

    def run_transeq(self):
        nul_seq = self.work_dir + "/" + self.option('sample') + ".CDS.fasta"
        port_seq = self.work_dir + "/" + self.option('sample') + ".faa"
        cmd = '{} -trim -table 11 -sequence {} -outseq {}'.format(self.transeq, nul_seq, port_seq)
        command = self.add_command("transeq", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("翻译蛋白序列运行完成")
        else:
            self.set_error("翻译蛋白序列运行出错!")

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
        if os.path.exists(self.work_dir + '/all.fasta'):
            os.remove(self.work_dir + '/all.fasta')
        cmd = "{} blastp -d {} -q {} -o {} --max-target-seqs 1".format(self.diamond, self.index,self.all_core_gene_fa,self.work_dir +'/' + self.option('sample') + ".matches.m8")
        self.logger.info(cmd)
        self.logger.info("开始运行run_get_gene")
        command = self.add_command("run_get_gene", cmd)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("运行run_get_gene完成")
        else:
            self.set_error("运行run_get_gene运行出错!")

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

    def run_stat(self):
        """
        合并16s和housekeeping结果
       :return:
        """
        hous_num = 0
        self.sample = self.option('sample')
        if os.path.exists(self.work_dir+"/"+self.sample + ".matches.m8"):
            with open(self.work_dir+"/"+self.sample + ".matches.m8", "r") as f:
                lines =f.readlines()
                hous_num =len(lines)
        s16_num = 0
        if os.path.exists(self.work_dir+"/"+self.sample + ".16s.fasta"):
            s16_num =len(list(SeqIO.parse(self.work_dir+"/"+self.sample + ".16s.fasta", "fasta")))
        with open(self.output_dir+"/"+self.sample + ".stat.xls", "w") as f:
            f.write("{}\t{}\t{}\n".format(self.sample, s16_num, hous_num))

    def set_output(self):
        if os.path.exists(self.output_dir+"/"+self.sample + ".16S.ffn"):
            os.remove(self.output_dir + "/" + self.sample + ".16S.ffn")
        if os.path.exists(self.work_dir + "/" + self.sample + ".16s.fasta"):
            os.link(self.work_dir + "/" + self.sample + ".16s.fasta", self.output_dir+ "/" + self.sample + ".16S.ffn")
            self.option("s16", self.output_dir + "/" + self.sample + ".16S.ffn")
        self.option("house", self.output_dir + "/" + self.sample + ".core_gene.fa")
        self.option("stat", self.output_dir + "/" + self.sample + ".stat.xls")

    def run(self):
        super(GffGenomePredictTool, self).run()
        self.run_getseq()
        self.run_transeq()
        self.run_index()
        self.run_get_gene()
        self.run_get_seq()
        self.run_stat()
        self.set_output()
        self.end()
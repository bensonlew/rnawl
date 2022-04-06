#!/usr/bin/env python
# -*- coding: utf-8 -*-
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os,re
from Bio import SeqIO
from mbio.packages.align.blast.xml2table  import   xml2table_coverage

class GetHgeneTreeAgent(Agent):
    """
    用于获取看家基因进化树序列
    version 1.0
    author: gaohao
    last_modify: 2018.04.21
    """

    def __init__(self, parent):
        super(GetHgeneTreeAgent, self).__init__(parent)
        options = [
            {"name": "gene_gff", "type": "infile", "format": "gene_structure.gff3"},  #
            {"name": "seq_faa", "type": "infile", "format": "sequence.fasta"},  #
            {"name": "seq_fa", "type": "infile", "format": "sequence.fasta"},  # gene核酸文件
            {"name": "out", "type": "outfile", "format": "sequence.fasta"},
            {"name": "sample_name", "type": "string"},
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
        self._cpu = 8
        self._memory = '10G'

    def end(self):
        super(GetHgeneTreeAgent, self).end()


class GetHgeneTreeTool(Tool):
    """
    version 1.0
    """
    def __init__(self, config):
        super(GetHgeneTreeTool, self).__init__(config)
        self.gene =self.option('seq_faa').prop['path']
        self.python = "/miniconda2/bin/python"
        self.perl_path = "/program/perl-5.24.0/bin/perl"
        self.perl_script = self.config.PACKAGE_DIR + "/bacgenome/"
        self.diamond = "/bioinfo/align/diamond-0.8.35/diamond"
        self.core_gene = self.config.SOFTWARE_DIR + "/database/bacgenome/House-keeping_gene/20210421/Core_gene"
        self.all_core_gene = self.config.SOFTWARE_DIR + "/database/bacgenome/House-keeping_gene/20210421/all_Coregene"
        self.all_core_gene_fa = self.config.SOFTWARE_DIR + "/database/bacgenome/House-keeping_gene/20210421/all.Coregene.name.faa"
        self.list = self.config.SOFTWARE_DIR + "/database/bacgenome/House-keeping_gene/20210421/cor_gene.list"

    def run_get_gene(self):
        if os.path.exists(self.work_dir + '/all.fasta'):        #guanqing.zou 20180904
            os.remove(self.work_dir + '/all.fasta')
        cmd = "{} blastp -d {} -q {} -o {} -p 8".format(self.diamond, self.core_gene,self.gene,self.work_dir +'/' + self.option('sample_name') + ".matches.m8")
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
            self.set_error("运行run_get_seq运行出错!", code="31402002")

    def run_diamond(self):
        cmd = "{} blastp -d {} -q {} -o {}  --max-target-seqs 19".format(self.diamond, self.all_core_gene, self.work_dir +  '/' + self.option('sample_name') + '.cor_gene.fa',
                                                   self.work_dir + '/' + self.option('sample_name') + ".last.m5")
        cmd += ' -f 5'  #zouguanqing 20190328
        cmd += ' -p 8'
        self.logger.info(cmd)
        self.logger.info("开始运行run_diamond")
        command = self.add_command("run_diamond", cmd)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("运行run_diamond完成")
        else:
            self.set_error("运行run_diamond运行出错!", code="31402003")
        xml2table_coverage(self.work_dir + '/' + self.option('sample_name') + ".last.m5", self.work_dir + '/' + self.option('sample_name') + ".last.m8") #zouguanqing 20190328

    def run_get_geneseq(self):
        cmd = "{} {}combin_fa.pl {} {} ".format(self.perl_path, self.perl_script,self.all_core_gene_fa , self.work_dir + '/' + self.option('sample_name') + ".last.m8")
        self.logger.info(cmd)
        self.logger.info("开始运行run_get_geneseq")
        command = self.add_command("run_get_geneseq", cmd)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            os.system('cat %s %s >%s' %(self.work_dir + '/' + self.option('sample_name') + '.cor_gene.fa',self.work_dir + '/all.fasta',self.work_dir + '/all.hgene.fa'))
            self.logger.info("运行run_get_geneseq完成")
        else:
            self.set_error("运行run_get_geneseq运行出错!", code="31402004")

    def get_balst_fungene(self):
        """
        主要是根据gene比对提取基因的详细信息
        :return:
        """
        cmd = "{} {}coregene_gff.py -b {} -g {} -o {}".format(self.python, self.perl_script, self.work_dir +'/' + self.option('sample_name') + ".matches.m8", self.option("gene_gff").prop['path'], self.work_dir+"/"+ self.option('sample_name') + ".coregene.m8")
        self.logger.info(cmd)
        self.logger.info("开始运行get_balst_fungene")
        command = self.add_command("get_balst_fungene", cmd)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("运行get_balst_fungene完成")
        else:
            self.set_error("运行get_balst_fungene运行出错!")

    def set_output(self):
        path =self.work_dir + '/all.hgene.fa'
        if os.path.exists(path):
            self.option('out').set_path(path)

    def get_fa(self):
        dict = {}
        with open(self.work_dir + '/' + self.option('sample_name') + ".matches.m8","r") as f:
            lines = f.readlines()
            for line in lines:
                lin = line.strip().split("\t")
                dict[lin[0]] = lin[1]
        self.getchuli(dict,self.option("seq_fa").prop['path'], self.work_dir + '/' + self.option('sample_name') + ".cor_gene.fnn")
        self.getchuli(dict, self.option("seq_faa").prop['path'], self.work_dir + '/' + self.option('sample_name') + ".cor_gene.faa")

    def getchuli(self,dict, file, out):
        list1 = []
        for fa_iterator in SeqIO.parse(file, "fasta"):
            id = fa_iterator.id
            if id in dict.keys():
                fa_iterator.id = id + "_" + dict[id]
                list1.append(fa_iterator)
        SeqIO.write(list1, out, "fasta")

    def run(self):
        """
        运行
        """
        super(GetHgeneTreeTool, self).run()
        self.run_get_gene()
        if os.path.getsize(self.work_dir +'/' + self.option('sample_name') + ".matches.m8") != 0:
            self.run_get_seq()
            self.run_diamond()
            self.run_get_geneseq()
            self.get_balst_fungene()
            self.get_fa()
            self.set_output()
        self.end()
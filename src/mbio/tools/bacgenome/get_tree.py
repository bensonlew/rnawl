#!/usr/bin/env python
# -*- coding: utf-8 -*-
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os,re


class GetTreeAgent(Agent):
    """
    用于获取16s进化树序列
    version 1.0
    author: gaohao
    last_modify: 2018.04.20
    """

    def __init__(self, parent):
        super(GetTreeAgent, self).__init__(parent)
        options = [
            {"name": "rrna_gff", "type": "infile", "format": "gene_structure.gff3"},  #
            {"name": "genome_fa", "type": "infile", "format": "sequence.fasta"},  #
            {"name": "out", "type": "outfile", "format": "sequence.fasta"},
            {"name": "sample_name", "type": "string"},
        ]
        self.add_option(options)
        self.list =[]


    def check_options(self):
        """
        检测参数是否正确
        """
        if not self.option('rrna_gff').is_set:
            raise OptionError("请设置基因组基因rRNAde gff3不存在！", code="31402101")
        if not self.option('genome_fa').is_set:
            raise OptionError("请设置基因组基因组序列不存在！", code="31402102")

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 5
        self._memory = '10G'

    def end(self):
        super(GetTreeAgent, self).end()


class GetTreeTool(Tool):
    """
    version 1.0
    """
    def __init__(self, config):
        super(GetTreeTool, self).__init__(config)
        self.rrna =self.option('rrna_gff').prop['path']
        self.seq = self.option('genome_fa').prop['path']
        self.perl_path = "/program/perl-5.24.0/bin/perl"
        self.perl_script = self.config.PACKAGE_DIR + "/bacgenome/"
        self.blastn = "/bioinfo/align/ncbi-blast-2.3.0+/bin/blastn"
        #self.database = self.config.SOFTWARE_DIR + "/database/bacgenome/16S/16s"
        self.database = self.config.SOFTWARE_DIR + "/database/bacgenome/16S/20210421/16s"
        self.fa =self.config.SOFTWARE_DIR + "/database/bacgenome/16S/20210421/16s.fasta"

    def run_get_16s(self):
        cmd = "{} {}tiqu_16s.pl {} {} {}".format(self.perl_path, self.perl_script,self.seq, self.rrna, self.option('sample_name'))
        self.logger.info(cmd)
        self.logger.info("开始运行run_get_16s")
        command = self.add_command("run_get_16s", cmd)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("运行run_get_16s完成")
        else:
            self.set_error("运行run_get_16s运行出错!", code="31402101")


    def run_blastn(self):
        cmd = "{} -db {} -query {} -out {} -outfmt 6 -max_target_seqs 200".format(self.blastn, self.database, self.work_dir + '/16s.fasta', self.work_dir + '/ori_blast_out.xls')
        self.logger.info(cmd)
        self.logger.info("开始运行run_blastn")
        command = self.add_command("run_blastn", cmd)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("运行run_blastn完成")
        else:
            self.set_error("运行run_blastn运行出错!", code="31402102")

    def run_get_seq(self):
        cmd = "{} {}combin_fa.pl {} {}".format(self.perl_path, self.perl_script,self.fa , self.work_dir + '/blast_out.xls')
        self.logger.info(cmd)
        self.logger.info("开始运行run_get_seq")
        command = self.add_command("run_get_seq", cmd)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            os.system('cat %s %s >%s' %(self.work_dir + '/16s.fasta',self.work_dir + '/all.fasta',self.work_dir + '/16s.tree.fa'))
            self.logger.info("运行run_get_seq完成")
        else:
            self.set_error("运行run_get_seq运行出错!", code="31402103")

    def set_output(self):
        path =self.work_dir + '/16s.tree.fa'
        if os.path.exists(path):
            self.option('out').set_path(path)

    def pick_species(self):
        with open(self.work_dir + '/ori_blast_out.xls') as f, open(self.work_dir + '/blast_out.xls','w') as fw:
            lines = f.readlines()
            dict = {}
            num =0
            for line in lines:
                lin = line.strip().split("\t")
                id =lin[1].split("_rRNA")[0]
                if id not in dict:
                    num +=1
                    dict[id] = [num,lin]
            rule =[1,2,3,4,5]
            if len(dict.keys())<=20:
                for i in range(6,21):
                    rule.append(i)
            else:
                num1 =(len(dict.keys())-5)/15
                while len(rule) < 20:
                    rule.append(rule[-1] + num1)
            self.logger.info(rule)
            for i in rule:
                for key in dict.keys():
                    if i == dict[key][0]:
                        fw.write("\t".join(dict[key][1]) + "\n")


    def run(self):
        """
        运行
        """
        super(GetTreeTool, self).run()
        self.run_get_16s()
        self.run_blastn()
        self.pick_species()
        self.run_get_seq()
        self.set_output()
        self.end()

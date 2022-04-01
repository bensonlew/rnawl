# -*- coding: utf-8 -*-
# __author__ = 'zouguanqing'
# version 1.0
# last_modify: 2017.12.27

import os
import re
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import subprocess


class ProdigalAgent(Agent):
    """
    prodigal 进行基因预测
    """

    def __init__(self, parent):
        super(ProdigalAgent, self).__init__(parent)
        options = [
            {"name": "input_genome", "type": "infile", "format": "sequence.fasta"},  # 组装拼接好的scaffold文件
            {"name": "gene_faa", "type": "outfile", "format": "sequence.fasta"},
            {"name": "gene_ffn", "type": "outfile", "format": "sequence.fasta"},
            {"name": "gene_gff", "type": "outfile", "format": "gene_structure.gff3"},
            {"name": "sample", "type": "string","default":"out"},
            {"name": "trans_table","type":"string", "defalut":"11"},  #翻译的系统
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("input_genome").is_set:
            raise OptionError("必须设置参数input_genome")

    def set_resource(self):
        self._cpu = 2
        self._memory = '2G'

    def end(self):
        super(ProdigalAgent, self).end()


class ProdigalTool(Tool):
    def __init__(self, config):
        super(ProdigalTool, self).__init__(config)
        self.genome_fasta = self.option("input_genome").prop['path']
        #self.sample_name = (self.genome_fasta).split('/')[-1].rsplit('.')[0]
        #self.perl_path = "/program/perl-5.24.0/bin/perl"
        #self.perl_script = self.config.PACKAGE_DIR + "/sequence/scripts/"
        self.prodigal_path = self.config.SOFTWARE_DIR +"/bioinfo/metaGenomic/Prodigal-2.6.3/prodigal"

    def run_prodigal(self):
        cmd = "{0} -i {2}  -o {1}.gff  -f gff -a {1}.faa -d {1}.ffn ".format(self.prodigal_path,self.option('sample'),self.genome_fasta)
        try:
            self.logger.info(cmd)
            subprocess.check_output(cmd, shell=True)
            self.logger.info("prodigal运行完成")
        except subprocess.CalledProcessError:
            self.set_error("prodigal运行出错")


    def set_output(self):
        if os.path.exists(self.output_dir+'/'+self.option('sample')+'.gff'):
            os.remove(self.output_dir+'/'+self.option('sample')+'.gff')
        os.link(self.work_dir+'/'+self.option('sample')+'.gff', self.output_dir+'/'+self.option('sample')+'.gff')
        self.option('gene_gff',  self.output_dir+'/'+self.option('sample')+'.gff')
        self.option('gene_ffn',  self.work_dir+'/'+self.option('sample')+'.ffn')
        self.option('gene_faa',  self.work_dir+'/'+self.option('sample')+'.faa')


    def run(self):
        super(ProdigalTool, self).run()
        self.run_prodigal()
        self.set_output()
        self.end()

# -*- coding: utf-8 -*-

import os,re
import shutil
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from mbio.packages.ref_rna.trans_step import step_count
from mbio.packages.fungi_genome.get_genemark_key import update_gm_key_pip


class GenePredictsModule(Module):
    """
    progigal小工具基因预测
    """

    def __init__(self, work_id):
        super(GenePredictsModule, self).__init__(work_id)
        option = [
            {"name": "input_genome", "type": "infile", "format": "sequence.fasta"},  # 组装拼接好的scaffold文件
            {"name": "genome_type", "type": "string", "default": ""},  # scaffold类型是细菌序列还是质粒序列("plas"), default": "" 是染色体
            {"name": "seq", "type": "outfile", "format": "sequence.fasta"},  # 预测基因DNA序列
            {"name": "prot_seq", "type": "outfile", "format": "sequence.fasta"},  # 预测基因prot序列
            {"name": "orf_prefix", "type": "string", "default": ""},  # ORF，自定义前缀
            {"name": "gene_statistics", "type": "outfile", "format": "sequence.profile_table"},  # 预测基因统计文件
            {"name": "gff", "type": "outfile", "format": "gene_structure.gff3"},  # 预测生成的GFF格式
            {"name": "software_list", "type":"string", "default":""},  #使用的软件，逗号分割
            {"name": "trna", "type": "infile", "format": "gene_structure.gff3"},   # 用来取舍预测到的基因，和ncrna重叠的基因去除
            {"name": "rrna", "type": "infile", "format": "gene_structure.gff3"},   # 用来取舍预测到的基因，和ncrna重叠的基因去除
            {"name": "trans_code", "type" : "string","default":"11"}
        ]

        self.add_option(option)
        self.prodigal = self.add_tool('predict.prodigal')
        self.collect = self.add_tool('tool_lab.gene_tidy')
        self.tools = []




    def check_options(self):
        """
        检查参数
        :return:
        """
        if not self.option("input_genome").is_set:
            raise OptionError("必须设置参数input_genome")
        self.genome_fasta = self.option('input_genome').prop['path']
        self.sample_name = (self.genome_fasta).split('/')[-1].rsplit('.')[0]

    def run_prodigal(self):
        opts = {
            'input_genome': self.option('input_genome'),
        }
        self.prodigal.set_options(opts)
        self.prodigal.run()

    def run_collect(self):
        opts ={
            'genome':self.option('input_genome'),
            'sample_name' : self.sample_name
            }
        opts['prodigal'] = self.prodigal.option('gene_gff')
        self.collect.set_options(opts)
        self.collect.run()

    def run(self):
        super(GenePredictsModule, self).run()
        self.prodigal.on('end',self.run_collect)
        self.collect.on('end',self.set_output)
        self.run_prodigal()

    def set_output(self, event):

        self.logger.info("设置结果目录")
        files = [self.sample_name + ".fnn",self.sample_name + ".faa", "gene_statistics.xls", self.sample_name+'.predict.gff']
        for f in files:
            if os.path.exists(self.output_dir + '/' + f):
                os.remove(self.output_dir + '/' + f)
        for f in files:
            os.link(self.collect.output_dir + '/' + f, self.output_dir + '/' + f)

        nul_path = self.output_dir + "/" + self.sample_name + ".fnn"
        port_path = self.output_dir + "/" + self.sample_name + ".faa"
        statistics_path = self.output_dir + "/gene_statistics.xls"
        gff_path = self.output_dir + "/" + self.sample_name + ".predict.gff"
        if os.path.getsize(nul_path):
            self.option('seq', nul_path)
            self.option('prot_seq', port_path)
            step_count(self.option('seq').prop['path'], self.collect.work_dir + "/fnn_stat.xls", 11, 100,
                       self.output_dir + "/length_distribute.txt")
        self.option('gene_statistics', statistics_path)
        self.option('gff', gff_path)
        self.end()


    def end(self):
        super(GenePredictsModule, self).end()

# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'
# version 1.0
# last_modify: 20190918

import os,re
import shutil
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from mbio.packages.ref_rna.trans_step import step_count


class GenePredictCdsModule(Module):
    """
    细基因组预测多种软件汇总
    """

    def __init__(self, work_id):
        super(GenePredictCdsModule, self).__init__(work_id)
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
            {"name": "trans_code", "type" : "string","default":"11"},
            {"name": "genome_name", "type" : "string"}, #基因组名称用于以区别组装序列的名称
        ]

        self.add_option(option)
        self.glimmer = self.add_tool('bac_comp_genome.glimmer_and_genemark')
        self.genemark = self.add_tool('bac_comp_genome.glimmer_and_genemark')
        self.prodigal = self.add_tool('predict.prodigal')
        self.collect = self.add_tool('bac_comp_genome.gene_tidy')
        self.software = []
        self.tools = []


    def check_options(self):
        """
        检查参数
        :return:
        """
        if not self.option("input_genome").is_set:
            raise OptionError("必须设置参数input_genome")
        if not self.option("genome_name"):
            raise OptionError("必须设置参数genome_name")
        if self.option('software_list') == "":
            if self.option("genome_type") == "plas":
                self.software.append('genemark')
            else:
                self.software.append('glimmer')
        else:
            self.software = self.option('software_list').split(',')
        self.genome_fasta = self.option('input_genome').prop['path']
        self.sample_name = self.option('genome_name')


    def run_glimmer(self):
        opts = {
            'input_genome': self.option('input_genome'),
            'software' : 'glimmer',
            'orf_prefix': self.option('orf_prefix'),
            'genome_name': self.option('genome_name')
        }
        self.glimmer.set_options(opts)
        self.glimmer.run()

    def run_genemark(self):
        opts = {
            'input_genome': self.option('input_genome'),
            'software' : 'genemark',
            'orf_prefix': self.option('orf_prefix'),
            'genome_name': self.option('genome_name')
        }
        self.genemark.set_options(opts)
        self.genemark.run()

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
        if 'glimmer' in self.software:
            opts['glimmer'] = self.glimmer.option('gff')
        if 'genemark' in self.software:
            opts['genemark'] = self.genemark.option('gff')
        if 'prodigal' in self.software:
            opts['prodigal'] = self.prodigal.option('gene_gff')
        if self.option('trna').is_set:
            opts['trna'] = self.option('trna')
        if self.option('rrna').is_set:
            opts['rrna'] = self.option('rrna')
        self.collect.set_options(opts)
        self.collect.run()


    def run(self):
        super(GenePredictCdsModule, self).run()
        self.tools = []
        for tool in self.software:
            if tool == 'glimmer':
                self.tools.append(self.glimmer)
            elif tool == 'genemark':
                self.tools.append(self.genemark)
            elif tool == 'prodigal':
                self.tools.append(self.prodigal)
        if len(self.tools) == 1:
            self.tools[0].on('end',self.run_collect)
        else:
            self.on_rely(self.tools,self.run_collect)
        self.collect.on('end',self.set_output)


        for tool in self.software:
            if tool == 'glimmer':
                self.run_glimmer()
            elif tool == 'genemark':
                self.run_genemark()
            elif tool == 'prodigal':
                self.run_prodigal()


    def set_output(self):
        """
        设置结果目录
        :param event:
        :return:
        """
        self.logger.info("设置结果目录")
        if len(self.tools) == 1:
            for tool in self.tools:
                for file in os.listdir(tool.output_dir):
                    file_path = os.path.join(tool.output_dir, file)
                    new_path = os.path.join(self.output_dir, file)
                    if os.path.exists(new_path):
                        os.remove(new_path)
                    os.link(file_path, new_path)
                    if re.search(r'_CDS.gff', file):
                        gff_path = os.path.join(self.output_dir, file)
                        self.option('gff', gff_path)
                    if re.search(r'.faa', file):
                        faa_path = os.path.join(self.output_dir, file)
                        self.option('prot_seq', faa_path)
                    if re.search(r'.fnn', file):
                        fnn_path = os.path.join(self.output_dir, file)
                        self.option('seq', fnn_path)
        else:
            files = [self.sample_name + ".ffn",self.sample_name + ".faa", "gene_statistics.xls", self.sample_name+'.predict.gff']
            for f in files:
                if os.path.exists(self.output_dir + '/' + f):
                    os.remove(self.output_dir + '/' + f)
            for f in files:
                os.link(self.collect.output_dir + '/' + f, self.output_dir + '/' + f)

            nul_path = self.output_dir + "/" + self.sample_name + ".ffn"
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
        super(GenePredictCdsModule, self).end()

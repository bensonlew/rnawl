# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
# version 1.0
# last_modify: 2021.01.18

from biocluster.module import Module
import os
import shutil
from biocluster.core.exceptions import OptionError
from mbio.packages.metagbin.common_function import link_dir
from Bio import SeqIO

class GenomeAniModule(Module):
    """
    单个基因组的16s比对后的物种注释和ani分析
    """
    def __init__(self, work_id):
        super(GenomeAniModule, self).__init__(work_id)
        options = [
            {"name": "type", "type": "string", "default":"GTDB"},  # GTDB or custom
            {'name': 'input', 'type': 'infile', "format": "sequence.fasta"},  ##基因组的16s序列
            {'name': 'house', 'type': 'infile', "format": "sequence.fasta"},  ##基因组的看家基因序列
            {"name": "genome", "type": "infile", "format": "sequence.fasta"},  ##样品的基因组文件
            {'name': 'custom_fa', 'type': 'infile', "format": "sequence.fasta"},  ##custom时的数据库序列
            {"name": "ref_s16", "type": "outfile", "format": "sequence.fasta"},
            {"name": "ref_genomes", "type": "infile", "format": "sequence.fasta_dir"},  ##自定义数据库的基因组文件
            {"name": "cus_table", "type": "infile", "format": "sequence.profile_table"},
            {"name": "s16_table", "type": "infile", "format": "sequence.profile_table"},
            {"name": "sample", "type": "string", "default": "out"},
        ]
        self.add_option(options)
        self.s16_align = self.add_tool('toolapps.rrna_align')
        self.get_genome = self.add_tool('toolapps.get_genome')
        self.list =[]

    def check_options(self):
        """
        检查参数
        :return:
        """
        if not self.option("genome").is_set:
            raise OptionError("必须设置参数genome")

    def get_fasta(self, input, output):
        """
        从样品16s的序列中取最长
        :return:
        """
        dict = {}
        for i in SeqIO.parse(input, "fasta"):
            dict[i] = len(i.seq)
        dict2 = sorted(dict.items(), key=lambda kv: (kv[1], kv[0]), reverse=True)
        dict2[0][0].id = self.option('sample')
        list = []
        list.append(dict2[0][0])
        SeqIO.write(list, output, "fasta")

    def run_16s_align(self):
        self.get_fasta(self.option("input").prop['path'], self.work_dir+"/16s.fasta")
        opts = ''
        if self.option('type') in ['GTDB']:
            opts = {
                'type': self.option('type'),
                'input': self.work_dir+"/16s.fasta",
                'sample': self.option('sample'),
            }
        elif self.option('type') in ['custom']:
            opts = {
                'type': self.option('type'),
                'input': self.work_dir+"/16s.fasta",
                'custom_fa': self.option("custom_fa"),
                'sample': self.option('sample'),
            }
        self.s16_align.set_options(opts)
        self.s16_align.run()

    def run_getgenome(self):
        if self.option("cus_table").is_set:
            opts = {
                'blast': self.s16_align.option('blast'),
                'cus_table': self.option('cus_table'),
            }
        else:
            opts = {
                'blast': self.s16_align.option('blast'),
            }
        self.get_genome.set_options(opts)
        self.get_genome.run()

    def run_rrna_tree(self):
        self.s16_tree = self.add_tool('toolapps.rrna_tree')
        self.list.append(self.s16_tree)
        os.system("cat {} {} >{}".format(self.get_genome.option("s16_fa").prop['path'], self.work_dir+"/16s.fasta", self.work_dir+"/all.16s.fasta"))
        opts = {
            'sample': self.option('sample'),
            's16_fa': self.work_dir+"/all.16s.fasta",
        }
        self.s16_tree.set_options(opts)
        self.s16_tree.run()

    def run_house_tree(self):
        self.house_tree = self.add_module('toolapps.house_tree')
        self.list.append(self.house_tree)
        opts = {
            'input': self.option("house"),
            'genomes': self.get_genome.option('genomes'),
            'sample': self.option('sample'),
        }
        self.house_tree.set_options(opts)
        self.house_tree.run()

    def run_genomeani(self):
        self.genome_ani = self.add_tool('toolapps.genome_ani')
        self.list.append(self.genome_ani)
        if self.option("cus_table").is_set:
            opts = {
                'sample': self.option('sample'),
                'genome': self.option('genome'),
                'blast': self.get_genome.option('blast_out'),
                'cus_table': self.option('cus_table'),
            }
        else:
            opts = {
                'sample': self.option('sample'),
                'genome': self.option('genome'),
                'blast': self.get_genome.option('blast_out'),
            }
        self.genome_ani.set_options(opts)
        self.genome_ani.run()

    def run_ani(self):
        if not self.get_genome.option('blast_out').is_set:
            self.run_rrna_tree()
            self.run_house_tree()
            self.on_rely(self.list, self.set_output)
        else:
            self.run_rrna_tree()
            self.run_house_tree()
            self.run_genomeani()
            self.on_rely(self.list, self.set_output)


    def run(self):
        """
        运行
        :return:
        """
        super(GenomeAniModule, self).run()
        self.s16_align.on("end", self.run_getgenome)
        self.get_genome.on("end", self.run_ani)
        self.run_16s_align()


    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        if not self.get_genome.option('blast_out').is_set:
            if os.path.exists(self.output_dir+"/" +self.option("sample")+".taxon.xls"):
                os.remove(self.output_dir+"/" +self.option("sample")+".taxon.xls")
            os.link(self.get_genome.output_dir+"/" +self.option("sample")+".taxon.xls", self.output_dir+"/" +self.option("sample")+".taxon.xls")
        else:
            if os.path.exists(self.output_dir+"/" +self.option("sample")+".taxon.xls"):
                os.remove(self.output_dir+"/" +self.option("sample")+".taxon.xls")
            os.link(self.genome_ani.output_dir+"/" +self.option("sample")+".taxon.xls", self.output_dir+"/" +self.option("sample")+".taxon.xls")
            if os.path.exists(self.output_dir+"/" +self.option("sample")+".anno.xls"):
                os.remove(self.output_dir+"/" +self.option("sample")+".anno.xls")
            os.link(self.genome_ani.output_dir+"/" +self.option("sample")+".anno.xls", self.output_dir+"/" +self.option("sample")+".anno.xls")
        if os.path.exists(self.output_dir + "/" + self.option("sample") + ".16s.nwk"):
            os.remove(self.output_dir + "/" + self.option("sample") + ".16s.nwk")
        os.link(self.s16_tree.output_dir + "/" + self.option("sample") + ".16s.nwk",
                self.output_dir + "/" + self.option("sample") + ".16s.nwk")
        if os.path.exists(self.output_dir + "/" + self.option("sample") + ".house_keeping.nwk"):
            os.remove(self.output_dir + "/" + self.option("sample") + ".house_keeping.nwk")
        os.link(self.house_tree.output_dir + "/" + self.option("sample") + ".house_keeping.nwk",
                self.output_dir + "/" + self.option("sample") + ".house_keeping.nwk")
        self.option("s16_table",self.s16_align.option("blast").prop['path'])
        self.end()

    def end(self):
        super(GenomeAniModule, self).end()
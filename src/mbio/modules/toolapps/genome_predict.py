# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
# version 1.0
# last_modify: 2021.01.05

from biocluster.module import Module
import os
import shutil
from biocluster.core.exceptions import OptionError
from Bio import SeqIO

class GenomePredictModule(Module):
    """
    单个基因组的16srRNA和housekeeping基因预测
    """
    def __init__(self, work_id):
        super(GenomePredictModule, self).__init__(work_id)
        options = [
            {"name": "input_genome", "type": "infile", "format": "sequence.fasta"},  # 组装拼接好的scaffold文件
            {"name": "s16", "type": "outfile", "format": "sequence.fasta"},
            {"name": "house", "type": "outfile", "format": "sequence.fasta"},
            {'name': 'stat', 'type': 'outfile', "format": "sequence.profile_table"},
        ]
        self.add_option(options)
        self.step.add_steps('rrna', 'housekeeping')
        self.rrna = self.add_tool('toolapps.genome_rrna')
        self.housekeeping = self.add_tool('toolapps.genome_housekeeping')
        self.list =[self.rrna, self.housekeeping]

    def check_options(self):
        """
        检查参数
        :return:
        """
        if not self.option("input_genome").is_set:
            raise OptionError("必须设置参数input_genome")

    def run_rrna(self):
        opts = {
            'input_genome': self.option('input_genome'),
        }
        self.rrna.set_options(opts)
        self.rrna.on('end', self.set_output, 'rrna')
        self.rrna.run()

    def run_housekeeping(self):
        self.sample = os.path.basename(self.option('input_genome').prop['path']).split(".fasta")[0]
        opts = {
            'input_genome': self.option('input_genome'),
            'sample': self.sample,
        }
        self.housekeeping.set_options(opts)
        self.housekeeping.on('end', self.set_output, 'housekeeping')
        self.housekeeping.run()

    def run_stat(self):
        """
        合并16s和housekeeping结果
       :return:
        """
        hous_num = 0
        if os.path.exists(self.housekeeping.output_dir+"/"+self.sample + ".matches.m8"):
            with open(self.housekeeping.output_dir+"/"+self.sample + ".matches.m8", "r") as f:
                lines =f.readlines()
                hous_num =len(lines)
        s16_num = 0
        if self.rrna.option("seq").is_set:
            s16_num =len(list(SeqIO.parse(self.rrna.option("seq").prop['path'], "fasta")))
        with open(self.output_dir+"/"+self.sample + ".stat.xls", "w") as f:
            f.write("{}\t{}\t{}\n".format(self.sample, s16_num, hous_num))
        if os.path.exists(self.output_dir+"/"+ self.sample + ".16S.ffn"):
            self.option("s16", self.output_dir+"/"+ self.sample + ".16S.ffn")
        self.option("house", self.output_dir + "/" + self.sample + ".core_gene.fa")
        self.option("stat", self.output_dir+"/"+self.sample + ".stat.xls")
        self.end()

    def run(self):
        """
        运行
        :return:
        """
        super(GenomePredictModule, self).run()
        self.on_rely(self.list, self.run_stat)
        self.run_housekeeping()
        self.run_rrna()


    def set_output(self, event):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        obj = event['bind_object']
        if event['data'] == 'rrna':
            if self.rrna.option("seq").is_set:
                if os.path.exists(self.output_dir+"/"+ self.sample + ".16S.ffn"):
                    os.remove(self.output_dir+"/"+ self.sample + ".16S.ffn")
                os.link(self.rrna.option("seq").prop['path'], self.output_dir+"/"+ self.sample + ".16S.ffn")
        elif event['data'] == 'housekeeping':
            if os.path.exists(self.output_dir + "/" + self.sample + ".core_gene.fa"):
                os.remove(self.output_dir + "/" + self.sample + ".core_gene.fa")
            os.link(self.housekeeping.output_dir + "/" + self.sample + ".core_gene.fa", self.output_dir + "/" + self.sample + ".core_gene.fa")

    def end(self):
        super(GenomePredictModule, self).end()
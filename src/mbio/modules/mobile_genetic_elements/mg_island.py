# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
import os
import shutil
import time
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from mbio.packages.metagbin.common_function import link_dir,link_file
from Bio import SeqIO

class MgIslandModule(Module):
    """
    单个基因组预测Island的预测，主要是三款软件的预测，将三款软件的预测结果整合，并处理
    author: gaohao
    last_modify: 2020.07.15
    """
    def __init__(self, work_id):
        super(MgIslandModule, self).__init__(work_id)
        options = [
            {"name": "genome_fa", "type": "infile", "format": "sequence.fasta"},  # 细菌的基因组fasta文件或宏基因组组装的fasta文件
            {"name": "gene_fna", "type": "infile",  "format": "sequence.fasta"},  # 预测的基因核酸文件
            {"name": "faa", "type": "infile", "format": "sequence.fasta"},  # 细菌的基因组fasta文件或宏基因组组装的fasta文件
            {"name": "ptt", "type": "string"},  # 类似gff
            {"name": "rnt", "type": "infile", "format": "metagbin.file_table"},  # rNA的统计文件
            {"name": "sample", "type": "string"},
        ]
        self.gihunter = self.add_module('mobile_genetic_elements.mg_gihunter_predict')
        self.dimob = self.add_module('mobile_genetic_elements.mg_island_dimob')
        self.islander = self.add_tool('mobile_genetic_elements.mg_island_islander')
        self.islander_stat = self.add_tool('mobile_genetic_elements.island_stat')
        self.add_option(options)

    def check_options(self):
        """
        检查参数
        :return:
        """
        pass

    def run_gihunter(self):
        opts = {
            "fna": self.option("genome_fa"),
            "ptt": self.option("ptt"),
            "rnt": self.option("rnt").prop['path'],
        }
        self.gihunter.set_options(opts)
        self.gihunter.run()

    def run_dimob(self):
        opts = {
            "sample_name":self.option("sample"),
            "gene_fna": self.option("gene_fna"),
            "faa": self.option("faa"),
            "ptt": self.option("ptt"),
            "genome_fa": self.option("genome_fa"),
        }
        self.dimob.set_options(opts)
        self.dimob.run()

    def run_islander(self):
        self.get_seq(self.option("genome_fa").prop['path'])
        opts = {
            "fa_dir": self.work_dir + "/seq_dir",
        }
        self.islander.set_options(opts)
        self.islander.run()

    def run_island_stat(self):
        if os.path.exists(self.work_dir + "/temp"):
            shutil.rmtree(self.work_dir + "/temp")
        os.mkdir(self.work_dir + "/temp")
        if self.dimob.option('out').is_set:
            if os.path.getsize(self.dimob.option('out').prop['path']) > 0:
                file = os.path.basename(self.dimob.option('out').prop['path'])
                os.link(self.dimob.option('out').prop['path'], self.work_dir + "/temp/" + file)
        if self.gihunter.option('out').is_set:
            file = os.path.basename(self.gihunter.option('out').prop['path'])
            os.link(self.gihunter.option('out').prop['path'], self.work_dir + "/temp/" + file)
        if self.islander.option('out').is_set:
            file = os.path.basename(self.islander.option('out').prop['path'])
            os.link(self.islander.option('out').prop['path'], self.work_dir + "/temp/" + file)
        if len(os.listdir(self.work_dir + "/temp")) >=1:
            opts = {
                "genome": self.option("genome_fa"),
                "dir": self.work_dir + "/temp",
                "sample_name": self.option("sample"),
                "gff": self.option("ptt")
            }
            self.islander_stat.set_options(opts)
            self.islander_stat.on("end", self.set_output)
            self.islander_stat.run()
        else:
            self.end()


    def run(self):
        """
        运行
        :return:
        """
        super(MgIslandModule, self).run()
        self.on_rely([self.islander, self.gihunter, self.dimob], self.run_island_stat)
        self.run_islander()
        self.run_dimob()
        self.run_gihunter()

    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        link_file(self.islander_stat.output_dir + "/" + self.option("sample") + ".GI_summary.xls", self.output_dir + "/" + self.option("sample") + ".island.xls")
        self.end()

    def get_seq(self, fa):
        if os.path.exists(self.work_dir + "/seq_dir"):
            shutil.rmtree(self.work_dir + "/seq_dir")
        os.mkdir(self.work_dir + "/seq_dir")
        for sercord in SeqIO.parse(fa, "fasta"):
            seq =self.work_dir + "/seq_dir/" + sercord.id + ".fna"
            SeqIO.write(sercord, seq, "fasta")

    def end(self):
        super(MgIslandModule, self).end()
# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
import os
import re
import time
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from mbio.packages.metagbin.common_function import link_dir,link_file

class MgIslandDimobModule(Module):
    """
    单个基因组预测island的预测，主要分为两部分，输入文件的处理和island的预测
    author: gaohao
    last_modify: 2020.07.15
    """
    def __init__(self, work_id):
        super(MgIslandDimobModule, self).__init__(work_id)
        options = [
            {"name": "gene_fna", "type": "infile", "format": "sequence.fasta"},  # 预测的基因核酸文件
            {"name": "faa", "type": "infile", "format": "sequence.fasta"},  # 细菌的基因组fasta文件或宏基因组组装的fasta文件
            {"name": "ptt", "type": "string"},  # 类似gff
            {"name": "genome_fa", "type": "infile", "format": "sequence.fasta"},
            {"name": "sample_name", "type": "string"},
            {"name": "out", "type": "outfile", "format": "sequence.profile_table"},
        ]
        self.island_dir = self.add_tool('mobile_genetic_elements.island_dir')
        self.island_dimob = self.add_tool('mobile_genetic_elements.mg_island_dimob')
        self.add_option(options)
        self.modules = []

    def check_options(self):
        """
        检查参数
        :return:
        """
        pass

    def run_file(self):
        opts = {
            "genome_fa":self.option("genome_fa"),
            "gene_fa": self.option("gene_fna"),
            "gene_faa": self.option("faa"),
            "gene_gff": self.option("ptt"),
        }
        self.island_dir.set_options(opts)
        self.island_dir.on("end", self.run_island)
        self.island_dir.run()

    def run_island(self):
        dir = self.island_dir.output_dir
        opts = {
            "sample_name": self.option("sample_name"),
            "dir": dir
        }
        self.island_dimob.set_options(opts)
        self.island_dimob.on("end", self.set_output)
        self.island_dimob.run()

    def run(self):
        """
        运行
        :return:
        """
        super(MgIslandDimobModule, self).run()
        self.run_file()

    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        if self.island_dimob.option('out').is_set:
            link_file(self.island_dimob.option('out').prop['path'], self.output_dir + "/" + self.option("sample_name") + ".Island_Dimob.xls")
            self.option("out", self.output_dir + "/" + self.option("sample_name") + ".Island_Dimob.xls")
        self.end()

    def end(self):
        super(MgIslandDimobModule, self).end()


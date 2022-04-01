# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
import os
import re
import time
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from mbio.packages.metagbin.common_function import link_dir

class IsPredictModule(Module):
    """
    单个基因组预测IS的程序，主要分为两部分，输入文件的处理和IS的预测
    author: gaohao
    last_modify: 2020.07.07
    """
    def __init__(self, work_id):
        super(IsPredictModule, self).__init__(work_id)
        options = [
            {"name": "genome_fa", "type": "infile", "format": "sequence.fasta"},  # 细菌的基因组fasta文件或宏基因组组装的fasta文件
            {"name": "gene_faa", "type": "infile",  "format": "sequence.fasta"},  # 预测的基因核酸文件
            {"name": "gene_fna", "type": "infile", "format": "sequence.fasta"},   # 预测的基因蛋白文件
            {"name": "gene_gff", "type": "string"},   # 预测的基因gff文件
        ]
        self.is_files = self.add_tool('mobile_genetic_elements.is_gene_dir')
        self.is_predict = self.add_tool('mobile_genetic_elements.is_predict')
        self.add_option(options)

    def run_file(self):
        opts = {
            "genome_fa": self.option("genome_fa"),
            "gene_faa": self.option("gene_faa"),
            "gene_fna": self.option("gene_fna"),
            "gene_gff": self.option("gene_gff"),
        }
        self.is_files.set_options(opts)
        self.is_files.on("end", self.run_is)
        self.is_files.run()

    def run_is(self):
        opts = {
            "input_fa":self.option("genome_fa"),
            "file_dir": self.is_files.output_dir,
            "gene_gff": self.option("gene_gff"),
        }
        self.is_predict.set_options(opts)
        self.is_predict.on("end", self.set_output)
        self.is_predict.run()

    def run(self):
        """
        运行
        :return:
        """
        super(IsPredictModule, self).run()
        self.run_file()

    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        if len(os.listdir(self.is_predict.output_dir)) > 0:
            link_dir(self.is_predict.output_dir, self.output_dir)
        self.end()

    def end(self):
        super(IsPredictModule, self).end()
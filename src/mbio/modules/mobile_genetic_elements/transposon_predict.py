# -*- coding: utf-8 -*-
# __author__ = 'gaohao'

import os
import shutil
import time
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from mbio.packages.metagbin.common_function import link_dir,link_file
from Bio import SeqIO

class TransposonPredictModule(Module):
    """
    可移动元件预测,主要是转座子和IS的数据兼容module
    author: gaohao
    last_modify: 2020.10.26
    """
    def __init__(self, work_id):
        super(TransposonPredictModule, self).__init__(work_id)
        options = [
            {"name": "genome_fa", "type": "infile", "format": "sequence.fasta"},  # 细菌的基因组fasta文件或宏基因组组装的fasta文件
            {"name": "gene_faa", "type": "infile", "format": "sequence.fasta"},  # 预测的基因核酸文件
            {"name": "gene_fna", "type": "infile", "format": "sequence.fasta"},  # 预测的基因蛋白文件
            {"name": "gene_gff", "type": "string"},  # 预测的基因gff文件
            {"name": "sample", "type": "string"},
        ]
        self.transposon = self.add_tool('mobile_genetic_elements.transposon_predict')
        self.is_pre = self.add_module('mobile_genetic_elements.is_predict')
        self.transposon_stat = self.add_tool('mobile_genetic_elements.transposon_combin')
        self.add_option(options)

    def check_options(self):
        """
        检查参数
        :return:
        """
        if not self.option("genome_fa").is_set:
            raise OptionError("必须设置参数基因组序列genome_fa")

    def run(self):
        """
        运行
        :return:
        """
        super(TransposonPredictModule, self).run()
        self.on_rely([self.transposon, self.is_pre], self.run_stat)
        self.run_transposon()
        self.run_is()


    def run_transposon(self):
        opts = {
            "input_fa": self.option("genome_fa"),
            "gene_gff": self.option("gene_gff"),
        }
        self.transposon.set_options(opts)
        self.transposon.run()

    def run_is(self):
        opts = {
            "genome_fa": self.option("genome_fa"),
            "gene_faa": self.option("gene_faa"),
            "gene_fna": self.option("gene_fna"),
            "gene_gff": self.option("gene_gff"),
        }
        self.is_pre.set_options(opts)
        self.is_pre.run()

    def run_stat(self):
        if os.path.exists(self.work_dir + "/temp"):
            shutil.rmtree(self.work_dir + "/temp")
        os.mkdir(self.work_dir + "/temp")
        if os.path.exists(self.transposon.output_dir + "/" + self.option("sample") + ".transposon.xls"):
            link_file(self.transposon.output_dir + "/" + self.option("sample") + ".transposon.xls", self.work_dir + "/temp/" + self.option("sample") + ".transposon.xls")
        if os.path.exists(self.is_pre.output_dir + "/" + self.option("sample") + ".is.xls"):
            link_file(self.is_pre.output_dir + "/" + self.option("sample") + ".is.xls",
                     self.work_dir + "/temp/" + self.option("sample") + ".is.xls")
        if len(os.listdir(self.work_dir + "/temp")) >0:
            opts = {
                "dir": self.work_dir + "/temp",
                "sample": self.option('sample')
            }
            self.transposon_stat.set_options(opts)
            self.transposon_stat.on("end", self.set_output)
            self.transposon_stat.run()
        else:
            self.end()


    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        if os.path.exists(self.output_dir + "/" + self.option("sample") + ".transposon.xls"):
            os.remove(self.output_dir + "/" + self.option("sample") + ".transposon.xls")
        if os.path.exists(self.transposon_stat.output_dir + "/" + self.option("sample") + ".transposon.xls"):
            os.link(self.transposon_stat.output_dir + "/" + self.option("sample") + ".transposon.xls", self.output_dir + "/" + self.option("sample") + ".transposon.xls")
        self.end()

    def end(self):
        super(TransposonPredictModule, self).end()
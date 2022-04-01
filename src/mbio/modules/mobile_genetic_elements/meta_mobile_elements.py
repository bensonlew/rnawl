# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
import os
import shutil
import time
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from mbio.packages.metagbin.common_function import link_dir,link_file
from Bio import SeqIO
import re

class MetaMobileElementsModule(Module):
    """
    宏基因组可移动元件的预测.分五部分，IS、island、integron、prompter、transposon
    author: gaohao
    last_modify: 2020.10.30
    """
    def __init__(self, work_id):
        super(MetaMobileElementsModule, self).__init__(work_id)
        options = [
            {"name": "genome_fa", "type": "infile", "format": "sequence.fasta"},  # 细菌的基因组fasta文件或宏基因组组装的fasta文件
            {"name": "gene_fna", "type": "infile",  "format": "sequence.fasta"},  # 预测的基因核酸文件
            {"name": "faa", "type": "infile", "format": "sequence.fasta"},  # 细菌的基因组fasta文件或宏基因组组装的fasta文件
            {"name": "ptt", "type": "string"},  # 类似gff
            {"name": "rnt", "type": "string"},  # rNA的统计文件
            {"name": "sample", "type": "string"},
        ]
        self.split_files = self.add_tool('mobile_genetic_elements.split_files')
        self.add_option(options)
        self.modules = []

    def check_options(self):
        """
        检查参数
        :return:
        """
        if not self.option("genome_fa").is_set:
            raise OptionError('必须输入genome_fa文件!')
        if not self.option("gene_fna").is_set:
            raise OptionError('必须输入gene_fna文件!')
        if not self.option("faa").is_set:
            raise OptionError('必须输入faa文件!')
        if not self.option("ptt"):
            raise OptionError('必须输入ptt文件!')
        if not self.option("rnt"):
            raise OptionError('必须输入rnt文件!')

    def run_split_files(self):
        opts = {
            "genome_fa": self.option("genome_fa"),
            "gene_faa": self.option("faa"),
            "gene_fna": self.option("gene_fna"),
            "gene_gff": self.option("ptt"),
            "gene_rnt": self.option("rnt"),
            "sample": self.option("sample"),
        }
        self.split_files.set_options(opts)
        self.split_files.run()

    def run_mobile(self):
        for dir in os.listdir(self.split_files.output_dir):
            mobile_elements = self.add_module('mobile_genetic_elements.meta_mobile')
            opts = {
                "genome_fa": self.split_files.output_dir + "/" + dir + "/" + dir + ".fna",
                "gene_fna": self.split_files.output_dir + "/" + dir + "/" + dir + ".ffn",
                "faa": self.split_files.output_dir + "/" + dir + "/" + dir + ".faa",
                "ptt": self.split_files.output_dir + "/" + dir + "/" + dir + ".ptt",
                "rnt": self.split_files.output_dir + "/" + dir + "/" + dir + ".rnt",
                "sample": dir
            }
            mobile_elements.set_options(opts)
            self.modules.append(mobile_elements)
        if len(self.modules) >1:
            self.on_rely(self.modules, self.set_output)
        elif len(self.modules) == 1:
            self.modules[0].on("end", self.set_output)
        for module in self.modules:
            module.run()

    def run(self):
        """
        运行
        :return:
        """
        super(MetaMobileElementsModule, self).run()
        self.split_files.on("end", self.run_mobile)
        self.run_split_files()

    def set_output(self, event):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        mge_file = []
        elem_file = []
        if len(self.modules) >1:
            for module in self.modules:
                files = os.listdir(module.output_dir)
                for file in files:
                    if re.search("mge.xls",file):
                        mge_file.append(module.output_dir + "/" + file)
                    elif re.search("element.xls",file):
                        elem_file.append(module.output_dir + "/" + file)
            if len(mge_file) >0:
                os.system("cat {} >{}".format(" ".join(mge_file), self.output_dir + "/" +self.option("sample") + "mge.xls"))
                os.system("cat {} >{}".format(" ".join(elem_file), self.output_dir + "/" +self.option("sample") + "element.xls"))
        elif len(self.modules) == 1:
            files = os.listdir(self.modules[0].output_dir)
            for file in files:
                if re.search("mge.xls", file):
                    link_file(self.modules[0].output_dir + "/" + file,
                              self.output_dir + "/" + self.option("sample") + "mge.xls")
                elif re.search("element.xls", file):
                    link_file(self.modules[0].output_dir + "/" + file,
                              self.output_dir + "/" + self.option("sample") + "element.xls")       
        self.end()

    def end(self):
        super(MetaMobileElementsModule, self).end()
# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'
import os
import re
import time
import shutil
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from mbio.packages.metagbin.common_function import link_dir

class IsPredictModule(Module):
    """
    单个基因组预测IS的程序，主要分为两部分，一是预测转座酶，二是is预测
    """
    def __init__(self, work_id):
        super(IsPredictModule, self).__init__(work_id)
        options = [
            {"name": "genome_fa", "type": "infile", "format": "sequence.fasta"},  # 细菌的基因组fasta文件或宏基因组组装的fasta文件
            {"name": "gene_faa", "type": "infile",  "format": "sequence.fasta"},  # 预测的基因蛋白文件
            {"name": "gene_fna", "type": "infile", "format": "sequence.fasta"},   # 预测的基因核酸文件
            {"name": "gene_gff", "type": "string"},   # 预测的基因gff文件
            {"name": "sample", "type": "string"}, ## 传入的样本名称
        ]
        self.is_files = self.add_tool('bacgenome.is_gene_dir')
        self.is_predict = self.add_tool('bacgenome.is_predict')
        self.enzyme = self.add_tool('bacgenome.diamond_enzyme')
        self.blast = self.add_tool("bacgenome.blast_is")
        self.add_option(options)

    def check_options(self):
        """
        检查参数
        :return:
        """
        pass

    def run_file(self):
        """
        format gene_faa、gene_gff
        """
        opts = {
            "genome_fa": self.option("genome_fa"),
            "gene_faa": self.option("gene_faa"),
            "gene_fna": self.option("gene_fna"),
            "gene_gff": self.option("gene_gff"),
            "sample": self.option("sample")
        }
        self.is_files.set_options(opts)
        self.is_files.on("end", self.run_is)
        self.is_files.on("end", self.run_pick_enzyme)
        self.is_files.run()

    def run_pick_enzyme(self):
        """
		运行diamond比对数据库提取enzyme的信息
        """
        opts = {
            "query": self.option("gene_faa"),
            "gff": os.path.join(self.is_files.output_dir, self.option("sample")+".gff"),
            "sample": self.option("sample"),
            "type": "transposase"
        }
        self.enzyme.set_options(opts)
        self.enzyme.run()


    def run_transposon(self):
        """
		预测转座子
        """
        opts = {
            "input_fa": self.option("genome_fa"),
            "gene_gff": os.path.join(self.is_files.output_dir, self.option("sample")+".gff"),
        }
        self.transposon.set_options(opts)
        self.transposon.run()

    def run_is(self):
        """
        运行ISEScan软件计算
        """
        opts = {
            "input_fa":self.option("genome_fa"),
            "file_dir": self.is_files.output_dir,
            "sample": self.option("sample")
        }
        self.is_predict.set_options(opts)
        self.is_predict.run()

    def run(self):
        """
        运行
        :return:
        """
        super(IsPredictModule, self).run()
        self.on_rely([self.enzyme, self.is_predict], self.run_blast_is)
        self.run_file()

    def run_blast_is(self):
        """
        运行blast软件比对数据库拿到is_name信息
        :return:
        """
        self.logger.info("开始用blast软件提取is的比对结果信息！")
        if os.path.exists(self.is_predict.output_dir + "/" + self.option("sample") + "_sequence.fna"):
            if os.path.getsize(self.is_predict.output_dir + "/" + self.option("sample") + "_sequence.fna") > 0:
                opts = {
                    "fna_file": self.is_predict.output_dir + "/" + self.option("sample") + "_sequence.fna",
                    "sample": self.option("sample")
                }
                self.blast.set_options(opts)
                self.blast.on("end", self.set_output)
                self.blast.run()
            else:
                self.set_output()
        else:
            self.set_output()

    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        sample_dir = os.path.join(self.output_dir, self.option("sample"))
        if os.path.exists(sample_dir):
            shutil.rmtree(sample_dir)
        os.mkdir(sample_dir)
        if len(os.listdir(self.is_predict.output_dir)) > 0:
            link_dir(self.is_predict.output_dir, sample_dir)
        if os.path.exists(sample_dir + "/" + self.option("sample") + ".enzyme.xls"):
            os.remove(sample_dir + "/" + self.option("sample") + ".enzyme.xls")
        if os.path.exists(self.enzyme.output_dir + "/" + self.option("sample") + ".enzyme.xls"):
            os.link(self.enzyme.output_dir + "/"+self.option("sample") + ".enzyme.xls", sample_dir + "/" + self.option("sample") + ".enzyme.xls")
        if os.path.exists(self.option("gene_fna").prop['path']):
            os.link(self.option("gene_fna").prop['path'], sample_dir + "/" + self.option("sample") + ".gene.ffn")
        if os.path.exists(self.blast.output_dir + "/" + self.option("sample") + ".blast.xls"):
            os.link(self.blast.output_dir + "/" + self.option("sample") + ".blast.xls", sample_dir + "/" + self.option("sample") + ".blast.xls")
        self.end()


    def end(self):
        super(IsPredictModule, self).end()
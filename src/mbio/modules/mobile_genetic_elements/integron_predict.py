# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
import os
import re
import time
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from mbio.packages.metagbin.common_function import link_dir

class IntegronPredictModule(Module):
    """
    单个基因组预测integron的预测，主要分为两部分，输入文件的处理和integron的预测
    author: gaohao
    last_modify: 2020.07.15
    """
    def __init__(self, work_id):
        super(IntegronPredictModule, self).__init__(work_id)
        options = [
            {"name": "genome_fa", "type": "infile", "format": "sequence.fasta"},  # 细菌的基因组fasta文件或宏基因组组装的fasta文件
            {"name": "gene_faa", "type": "infile",  "format": "sequence.fasta"},  # 预测的基因核酸文件
            {"name": "gene_gff", "type": "string"},   # 预测的基因gff文件
        ]
        self.integron_files = self.add_tool('mobile_genetic_elements.integron_dir')
        self.integron_predict = self.add_tool('mobile_genetic_elements.integron_predict')
        self.modules = []
        self.add_option(options)

    def check_options(self):
        """
        检查参数
        :return:
        """
        pass

    def run_file(self):
        opts = {
            "genome_fa": self.option("genome_fa"),
            "gene_faa": self.option("gene_faa"),
            "gene_gff": self.option("gene_gff"),
        }
        self.integron_files.set_options(opts)
        self.integron_files.on("end", self.run_integron)
        self.integron_files.run()

    def run_integron(self):
        opts = {
            "input_fa": self.option("genome_fa"),
            "file_dir": self.integron_files.output_dir
        }
        self.integron_predict.set_options(opts)
        self.integron_predict.on("end", self.set_output)
        self.integron_predict.run()

    def run(self):
        """
        运行
        :return:
        """
        super(IntegronPredictModule, self).run()
        self.run_file()

    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        for i in os.listdir(self.integron_predict.output_dir):
            if os.path.exists(self.output_dir + "/" + i):
                os.remove(self.output_dir + "/" + i)
            os.link(self.integron_predict.output_dir + "/" + i, self.output_dir + "/" + i)
        self.end()

    def end(self):
        super(IntegronPredictModule, self).end()
#-*- coding: utf-8 -*-
import os
import re
import time
import shutil
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from mbio.packages.bacgenome.common import link_dir

class IntegronPredictModule(Module):
    """
    单个基因组预测integron的预测，主要分为两部分，输入文件的处理和integron的预测
    """
    def __init__(self, work_id):
        super(IntegronPredictModule, self).__init__(work_id)
        options = [
            {"name": "genome_fa", "type": "infile", "format": "sequence.fasta"},  # 细菌的基因组fasta文件或宏基因组组装的fasta文件
            {"name": "gene_faa", "type": "infile",  "format": "sequence.fasta"},  # 预测的基因核酸文件
            {"name": "gene_gff", "type": "string"},   # 预测的基因gff文件
            {"name": "sample", "type": "string"}, ## 传入的样本名称
        ]
        self.integron_files = self.add_tool('bacgenome.format_dir')
        self.integron_predict = self.add_tool('bacgenome.integron_predict')
        self.add_option(options)

    def check_options(self):
        """
        检查参数
        :return:
        """
        if not self.option("genome_fa").is_set:
            raise OptionError("必须输入基因组序列！")

    def run_file(self):
        """
        format文件
        :return:
        """
        opts = {
            "genome_fa": self.option("genome_fa"),
            "gene_faa": self.option("gene_faa"),
            "gene_gff": self.option("gene_gff"),
        }
        self.integron_files.set_options(opts)
        self.integron_files.on("end", self.run_integron)
        self.integron_files.run()

    def run_integron(self):
        """
        运行软件integron_Finder进行预测
        :return:
        """
        opts = {
            "input_fa": self.option("genome_fa"),
            "file_dir": self.integron_files.output_dir,
            "sample": self.option("sample")
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
        sample_dir = os.path.join(self.output_dir, self.option("sample"))
        if os.path.exists(sample_dir):
            shutil.rmtree(sample_dir)
            os.mkdir(sample_dir)
        else:
            os.mkdir(sample_dir)
        for i in os.listdir(self.integron_predict.output_dir):
            if os.path.exists(sample_dir + "/" + i):
                os.remove(sample_dir + "/" + i)
            os.link(self.integron_predict.output_dir + "/" + i, sample_dir + "/" + i)
        self.end()

    def end(self):
        super(IntegronPredictModule, self).end()
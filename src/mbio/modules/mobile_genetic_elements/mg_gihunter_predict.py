# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
import os
import re
import time
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from mbio.packages.metagbin.common_function import link_dir,link_file

class MgGihunterPredictModule(Module):
    """
    单个基因组预测island的预测，主要Gihunter软件的预测
    author: gaohao
    last_modify: 2020.09.01
    """
    def __init__(self, work_id):
        super(MgGihunterPredictModule, self).__init__(work_id)
        options = [
            {"name": "fna", "type": "infile", "format": "sequence.fasta"},  # 细菌的基因组fasta文件或宏基因组组装的fasta文件
            {"name": "ptt", "type": "string"},  # 类似gff
            {"name": "rnt", "type": "string"},  # rNA的统计文件
            {"name": "out", "type": "outfile", "format": "sequence.profile_table"},
        ]
        self.add_option(options)
        self.gihunter_dir = self.add_tool('mobile_genetic_elements.gihunter_dir')
        self.gihunter = self.add_tool('mobile_genetic_elements.mg_gihunter_predict')
        self.modules =[]

    def check_options(self):
        """
        检查参数
        :return:
        """
        pass

    def get_dir(self):
        opts = {
            "genome_fa": self.option("fna"),
            "rnt": self.option("rnt"),
            "gff": self.option("ptt")
        }
        self.gihunter_dir.set_options(opts)
        self.gihunter_dir.on("end", self.run_island)
        self.gihunter_dir.run()

    def run_island(self):
        opts = {
            "dir": self.gihunter_dir.output_dir,
        }
        self.gihunter.set_options(opts)
        self.gihunter.on("end", self.set_output)
        self.gihunter.run()

    def run(self):
        """
        运行
        :return:
        """
        super(MgGihunterPredictModule, self).run()
        self.get_dir()

    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        if self.gihunter.option("out").is_set:
            link_file(self.gihunter.option("out").prop['path'], self.output_dir + "/gihunter.xls")
            self.option("out", self.output_dir + "/gihunter.xls")
        self.end()

    def end(self):
        super(MgGihunterPredictModule, self).end()
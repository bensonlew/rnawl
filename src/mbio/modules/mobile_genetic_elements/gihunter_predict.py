# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
import os
import re
import time
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from mbio.packages.metagbin.common_function import link_dir,link_file

class GihunterPredictModule(Module):
    """
    单个基因组预测island的预测，主要Gihunter软件的预测
    author: gaohao
    last_modify: 2020.09.01
    """
    def __init__(self, work_id):
        super(GihunterPredictModule, self).__init__(work_id)
        options = [
            {"name": "fna", "type": "infile", "format": "sequence.fasta"},  # 细菌的基因组fasta文件或宏基因组组装的fasta文件
            {"name": "ptt", "type": "string"},  # 类似gff
            {"name": "rnt", "type": "string"},  # rNA的统计文件
            {"name": "out", "type": "outfile", "format": "sequence.profile_table"},
        ]
        self.add_option(options)
        self.gihunter_dir = self.add_tool('mobile_genetic_elements.gihunter_dir')
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
        for file in os.listdir(self.gihunter_dir.output_dir):
            self.gihunter = self.add_tool('mobile_genetic_elements.gihunter_predict')
            opts = {
                "fna": self.gihunter_dir.output_dir + "/" + file + "/" + file + ".fna",
                "rnt": self.gihunter_dir.output_dir + "/" + file + "/" + file + ".rnt",
                "ptt": self.gihunter_dir.output_dir + "/" + file + "/" + file + ".ptt",
            }
            self.gihunter.set_options(opts)
            self.modules.append(self.gihunter)
        if len(self.modules) > 1:
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
        super(GihunterPredictModule, self).run()
        self.get_dir()

    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        if os.path.exists(self.output_dir + "/gihunter.xls"):
            os.remove(self.output_dir + "/gihunter.xls")
        if len(self.modules) > 1:
            with open(self.work_dir + "/gihunter.xls", "w") as g:
                for module in self.modules:
                    if module.option('out').is_set:
                        if self.get_num(module.option('out').prop['path']) >= 1:
                            with open(module.option('out').prop['path'], "r") as f:
                                lines = f.readlines()
                                for line in lines:
                                    lin = line.strip().split("\t")
                                    g.write("{}\t{}\t{}\t{}\n".format(lin[0], lin[1], lin[2], lin[3]))
        elif len(self.modules) == 1:
            with open (self.work_dir + "/gihunter.xls", "w") as g:
                if self.modules[0].option('out').is_set:
                    if self.get_num(self.modules[0].option('out').prop['path']) >= 1:
                        with open(self.modules[0].option('out').prop['path'], "r") as f:
                            lines = f.readlines()
                            for line in lines:
                                lin = line.strip().split("\t")
                                g.write("{}\t{}\t{}\t{}\n".format(lin[0], lin[1], lin[2], lin[3]))
        link_file(self.work_dir + "/gihunter.xls", self.output_dir + "/gihunter.xls")
        if os.path.getsize(self.output_dir + "/gihunter.xls") >0:
            self.option("out", self.output_dir + "/gihunter.xls")
        self.end()

    def get_num(self, file):
        with open(file, "r") as f:
            lines =f.readlines()
            num = len(lines)
        return num

    def end(self):
        super(GihunterPredictModule, self).end()
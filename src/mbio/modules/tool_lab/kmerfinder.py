# -*- coding: utf-8 -*-
# __author__ = 'hao.gao'
# last_modify: 2021.04.29

import os,shutil,re
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from mbio.packages.bac_comp_genome.common_function import link_dir,link_file
import pandas as pd

class KmerfinderModule(Module):
    """
    多个样品的Kmerfinder分析
    """

    def __init__(self, work_id):
        super(KmerfinderModule, self).__init__(work_id)
        option = [
            {"name": "fasta_dir", "type": "infile", "format": "tool_lab.fasta_dir"},  # 基因组文件夹
            {"name": "species", "type": "string", },  # 物种名称
        ]
        self.add_option(option)
        self.all_modules = []

    def check_options(self):
        """
        检查参数
        :return:
        """
        if not self.option("fasta_dir").is_set:
            raise OptionError("必须提fasta_dir！")

        if not self.option("species"):
            raise OptionError("必须提供species！")
        self.samples, self.fastas = self.get_sample(self.option("fasta_dir").prop['path']+"/list.txt")
        return True

    def run_kmerfinder(self):
        for sample in self.samples:
            kmerfinder = self.add_tool('tool_lab.kmerfinder')
            opts = {
                "fasta": self.option("fasta_dir").prop['path'] + "/" + self.fastas[sample],
                "species": self.option("species"),
                "sample_name": sample,
            }
            kmerfinder.set_options(opts)
            self.all_modules.append(kmerfinder)
        if len(self.all_modules) > 1:
            self.on_rely(self.all_modules, self.set_output)
        elif len(self.all_modules) == 1:
            self.all_modules[0].on("end", self.set_output)
        for module in self.all_modules:
            module.run()

    def set_output(self):
        if len(self.all_modules) > 1:
            for module in self.all_modules:
                if os.listdir(module.output_dir) > 0:
                    for i in os.listdir(module.output_dir):
                        link_file(module.output_dir + "/" + i, self.output_dir + "/" + i)
        elif len(self.all_modules) == 1:
            if os.listdir(self.all_modules[0].output_dir) >0:
                for i in os.listdir(self.all_modules[0].output_dir):
                    link_file(self.all_modules[0].output_dir+"/"+i, self.output_dir+"/"+i)
        self.end()

    def get_sample(self,file):
        list = []
        dict = {}
        with open (file, "r") as f:
            lines = f.readlines()
            for line in lines[1:]:
                lin = line.strip().split("\t")
                list.append(lin[0])
                dict[lin[0]] = lin[1]
        return list,dict

    def run(self):
        super(KmerfinderModule, self).run()
        self.run_kmerfinder()

    def end(self):
        super(KmerfinderModule, self).end()
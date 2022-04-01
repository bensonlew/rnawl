# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
import os
import re
import time
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from mbio.packages.metagbin.common_function import link_dir,link_file

class IslandIslanderModule(Module):
    """
    单个基因组预测island的预测，主要island的预测
    author: gaohao
    last_modify: 2020.07.15
    """
    def __init__(self, work_id):
        super(IslandIslanderModule, self).__init__(work_id)
        options = [
            {"name": "fa_dir", "type": "infile", "format": "sequence.fasta_dir"},  # 预测的基因核酸文件
            {"name": "out", "type": "outfile", "format": "sequence.profile_table"},
        ]
        self.add_option(options)
        self.modules =[]

    def check_options(self):
        """
        检查参数
        :return:
        """
        pass

    def run_island(self):
        for file in os.listdir(self.option("fa_dir").prop['path']):
            self.island = self.add_tool('mobile_genetic_elements.island_islander')
            opts = {
                "fna": self.option("fa_dir").prop['path'] + "/" + file
            }
            self.island.set_options(opts)
            self.modules.append(self.island)
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
        super(IslandIslanderModule, self).run()
        self.run_island()

    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        if len(self.modules) > 1:
            with open(self.work_dir + "/islander.xls", "w") as g:
                for module in self.modules:
                    if module.option('out').is_set:
                        if self.get_num(module.option('out').prop['path']) >= 2:
                            with open(module.option('out').prop['path'], "r") as f:
                                lines = f.readlines()
                                for line in lines[1:]:
                                    lin = line.strip().split("\t")
                                    g.write("{}\t{}\t{}\t{}\n".format(lin[0], lin[1], lin[3], lin[4]))

        elif len(self.modules) == 1:
            with open(self.work_dir + "/islander.xls", "w") as g:
                if self.modules[0].option('out').is_set:
                    if self.get_num(self.modules[0].option('out').prop['path']) >= 2:
                        with open(self.modules[0].option('out').prop['path'], "r") as f:
                            lines = f.readlines()
                            for line in lines[1:]:
                                lin = line.strip().split("\t")
                                g.write("{}\t{}\t{}\t{}\n".format(lin[0], lin[1], lin[3], lin[4]))
        link_file(self.work_dir + "/islander.xls", self.output_dir + "/islander.xls")
        if os.path.getsize(self.output_dir + "/islander.xls") >0:
            self.option("out", self.output_dir + "/islander.xls")
        self.end()

    def get_num(self, file):
        with open(file, "r") as f:
            lines =f.readlines()
            num = len(lines)
        return num

    def end(self):
        super(IslandIslanderModule, self).end()
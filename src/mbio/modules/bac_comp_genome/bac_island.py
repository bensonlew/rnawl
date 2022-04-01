# -*- coding: utf-8 -*-
# __author__ = 'hao.gao'
# last_modify: 2019.10.09

import os,shutil
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from mbio.packages.bac_comp_genome.common_function import link_dir,link_file

class BacIslandModule(Module):
    """
    多个样品的基因组岛预测
    """

    def __init__(self, work_id):
        super(BacIslandModule, self).__init__(work_id)
        option = [
            {"name": "gff_path", "type": "string"},  # 合并CDS/tRNA/rRNA的gff路径
            {"name": "gene_path", "type": "string"},  # string类型的基因预测dir路径
            {"name": "genome_path", "type": "string"},  # string类型的组装序列dir路径
            {"name": "gbk_path", "type": "string"},  # string类型gbk文件dir路径
            {"name": "seq_path", "type": "string"},  # string类型的单个样品每条序列dir路径
            {"name": "sample_list", "type": "infile", "format": "sequence.profile_table"},#文件中有样品信息
        ]
        self.add_option(option)
        self.blastclust =self.add_tool("bac_comp_genome.blastclust")
        self.all_modules = []

    def check_options(self):
        """
        检查参数
        :return:
        """
        if not self.option("sample_list").is_set:
            raise OptionError("必须提供样品的信息文件！")
        if not self.option("gene_path"):
            raise OptionError("必须提供gene_path的路径！")
        if not self.option("genome_path"):
            raise OptionError("必须提供genome_path的路径！")
        if not self.option("gbk_path"):
            raise OptionError("必须提供gbk_path的路径！")
        if not self.option("seq_path"):
            raise OptionError("必须提供seq_path的路径！")
        return True

    def run_all_island(self):
        self.samples = self.get_sample(self.option("sample_list").prop['path'])
        for sample in self.samples:
            island = self.add_module('bac_comp_genome.island')
            opts = {
                "fa_dir": self.option("seq_path") + sample + "/",
                "gbk_dir": self.option("gbk_path") + sample + "/",
                "gff": self.option("gff_path") + "/" + sample + ".gff",
                "genome": self.option("genome_path") + sample + ".fna",
                "sample_name": sample,
            }
            island.set_options(opts)
            self.all_modules.append(island)
        if len(self.all_modules) > 1:
            self.on_rely(self.all_modules, self.run_blastclust)
        elif len(self.all_modules) == 1:
            self.all_modules[0].on("end", self.run_blastclust)
        for module in self.all_modules:
            module.run()

    def run_blastclust(self):
        if not os.path.exists(self.work_dir + "/tmp_island"):
            os.mkdir(self.work_dir + "/tmp_island")
        if len(self.all_modules) > 1:
            for module in self.all_modules:
                if len(os.listdir(module.output_dir)) >= 1:
                    link_dir(module.output_dir, self.work_dir + "/tmp_island")
        elif len(self.all_modules) == 1:
            if len(os.listdir(self.all_modules[0].output_dir)) >= 1:
                link_dir(self.all_modules[0].output_dir, self.work_dir + "/tmp_island")
        if len(os.listdir(self.work_dir + "/tmp_island")) >=1:
            self.blastclust.set_options({
                "dir": self.work_dir + "/tmp_island",
                "method": "island",
                "sample_list": self.option("sample_list")
            })
            self.blastclust.on("end", self.set_output)
            self.blastclust.run()
        else:
            self.end()

    def set_output(self):
        if len(os.listdir(self.work_dir + "/tmp_island")) >= 1:
            for module in self.all_modules:
                if len(os.listdir(module.output_dir)) >= 1:
                    link_dir(module.output_dir, self.output_dir)
            link_dir(self.blastclust.output_dir, self.output_dir)
        self.end()

    def get_sample(self,file):
        list = []
        with open (file, "r") as f:
            lines = f.readlines()
            for line in lines[1:]:
                lin = line.strip().split("\t")
                list.append(lin[0])
        return list

    def run(self):
        super(BacIslandModule, self).run()
        self.run_all_island()

    def end(self):
        super(BacIslandModule, self).end()
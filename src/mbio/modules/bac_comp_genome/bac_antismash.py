# -*- coding: utf-8 -*-
# __author__ = 'hao.gao'
# last_modify: 2019.09.24

import os,shutil,re
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from mbio.packages.bac_comp_genome.common_function import link_dir,link_file

class BacAntismashModule(Module):
    """
    多个样品的注释模块开发
    """

    def __init__(self, work_id):
        super(BacAntismashModule, self).__init__(work_id)
        option = [
            {"name": "genome_path", "type": "string"},  # string类型的基因组dir路径
            {"name": "gbk_path", "type": "string"},  # string类型的gbk文件dir路径
            {"name": "sample_list", "type": "infile", "format": "sequence.profile_table"},#文件中有样品信息
        ]
        self.add_option(option)
        self.all_modules = []

    def check_options(self):
        """
        检查参数
        :return:
        """
        if not self.option("sample_list").is_set:
            raise OptionError("必须提供样品的信息文件！")
        if not self.option("genome_path"):
            raise OptionError("必须提供genome_path的路径！")
        if not self.option("gbk_path"):
            raise OptionError("必须提供gbk_path的路径！")
        return True

    def run_all_antismash(self):
        self.des = self.get_sample(self.option("sample_list").prop['path'])
        for des in self.des:
            sample, type_ann =des.split("\t")
            anno = self.add_module('bac_comp_genome.antismash')
            opts = {
                "genome": self.option("genome_path") + sample + ".fna",
                "gbk_dir": self.option("gbk_path") + sample + "/",
                "sample_name": sample,
            }
            anno.set_options(opts)
            self.all_modules.append(anno)
        if len(self.all_modules) > 1:
            self.on_rely(self.all_modules, self.antismash_stat)
        elif len(self.all_modules) == 1:
            self.all_modules[0].on("end", self.antismash_stat)
        for module in self.all_modules:
            module.run()

    def antismash_stat(self):
        if os.path.exists(self.work_dir + "/tmp_antismash"):
            shutil.rmtree(self.work_dir + "/tmp_antismash")
        os.mkdir(self.work_dir + "/tmp_antismash")
        if len(self.all_modules) >1:
            for module in self.all_modules:
                if len(os.listdir(module.output_dir)) >= 1:
                    link_dir(module.output_dir, self.work_dir + "/tmp_antismash")
        elif len(self.all_modules) == 1:
            if len(os.listdir(self.all_modules[0].output_dir)) >= 1:
                link_dir(self.all_modules[0].output_dir, self.work_dir + "/tmp_antismash")
        self.n = self.get_num(self.work_dir + "/tmp_antismash")
        if self.n >=1:
            self.antismash_sum = self.add_tool("bac_comp_genome.antismash_sum")
            self.antismash_sum.set_options({
                "dir": self.work_dir + "/tmp_antismash",
                "sample_list": self.option("sample_list")
            })
            self.antismash_sum.on("end", self.set_output)
            self.antismash_sum.run()
        else:
            self.end()

    def set_output(self):
        for module in self.all_modules:
            if len(os.listdir(module.output_dir)) >=1:
                link_dir(module.output_dir, self.output_dir)
        link_dir(self.antismash_sum.output_dir, self.output_dir)
        self.end()

    def get_num(self, dir):
        files = os.listdir(dir)
        list =[]
        for file in files:
            if re.search(r".antismash_anno.xls", file):
                list.append(file)
        return len(list)

    def get_sample(self,file):
        list = []
        with open(file, "r") as f:
            lines = f.readlines()
            for line in lines[1:]:
                lin = line.strip().split("\t")
                list.append(lin[0] + "\t" + lin[1])
        return list

    def run(self):
        super(BacAntismashModule, self).run()
        self.run_all_antismash()

    def end(self):
        super(BacAntismashModule, self).end()
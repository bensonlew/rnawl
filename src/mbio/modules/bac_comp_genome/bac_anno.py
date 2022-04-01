# -*- coding: utf-8 -*-
# __author__ = 'hao.gao'
# last_modify: 2019.09.24

import os,shutil
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from mbio.packages.bac_comp_genome.common_function import link_dir,link_file,anno_kegg

class BacAnnoModule(Module):
    """
    多个样品的注释模块开发
    """

    def __init__(self, work_id):
        super(BacAnnoModule, self).__init__(work_id)
        option = [
            {"name": "path", "type": "string"},  # string类型的dir路径
            {"name": "gff_path", "type": "string"},  # 合并CDS/tRNA/rRNA的gff路径
            {"name": "sample_list", "type": "infile", "format": "sequence.profile_table"},#文件中有样品信息
        ]
        self.add_option(option)
        self.compare_anno =self.add_tool("bac_comp_genome.compare_anno")
        self.anno_modules = []

    def check_options(self):
        """
        检查参数
        :return:
        """
        if not self.option("sample_list").is_set:
            raise OptionError("必须提供样品的信息文件！")
        return True

    def run_all_anno(self):
        self.samples = self.get_sample(self.option("sample_list").prop['path'])
        for sample in self.samples:
            anno = self.add_module('bac_comp_genome.annotation')
            opts = {
                "gene_seq": self.option("path") + sample + "/" + sample + "_CDS.faa",
                "gff": self.option("gff_path") + "/" + sample + ".gff",
                "sample": sample,
            }
            anno.set_options(opts)
            self.anno_modules.append(anno)
        if len(self.anno_modules) > 1:
            self.on_rely(self.anno_modules, self.run_comp_anno)
        elif len(self.anno_modules) == 1:
            self.anno_modules[0].on("end", self.set_output)
        for module in self.anno_modules:
            module.run()

    def run_comp_anno(self):
        if not os.path.exists(self.work_dir + "/anno"):
            os.mkdir(self.work_dir + "/anno")
        for module in self.anno_modules:
            link_dir(module.output_dir, self.work_dir + "/anno")
        for sample in self.samples:
            os.renames(
                self.work_dir + "/anno/" + sample + "/" + sample + "_Gram+_SignalP.txt",
                self.work_dir + "/anno/" + sample + "/" + sample + "_Gram--_SignalP.txt")
        self.compare_anno.set_options({
            "dir": self.work_dir + "/anno",
        })
        self.compare_anno.on("end", self.run_map)
        self.compare_anno.run()

    def run_map(self):
        if os.path.exists(self.work_dir + "/all_kegg.xls"):
            os.remove(self.work_dir + "/all_kegg.xls")
        list = []
        for module in self.anno_modules:
            list.append(module.option("kegg").prop['path'])
        anno_kegg(list,self.work_dir + "/all_kegg.xls")
        self.kegg_map = self.add_tool("bac_comp_genome.kegg_graph_info")
        self.kegg_map.set_options({
           "annotable": self.work_dir + "/all_kegg.xls"
        })
        self.kegg_map.on("end", self.set_output)
        self.kegg_map.run()

    def set_output(self):
        for module in self.anno_modules:
            link_dir(module.output_dir, self.output_dir)
        link_dir(self.compare_anno.output_dir, self.output_dir + "/all_anno")
        link_file(self.kegg_map.output_dir + "/kegg_graph_info.xls", self.output_dir + "/all_anno/kegg/kegg_graph_info.xls")
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
        super(BacAnnoModule, self).run()
        self.run_all_anno()

    def end(self):
        super(BacAnnoModule, self).end()
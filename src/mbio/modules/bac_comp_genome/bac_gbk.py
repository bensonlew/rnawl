# -*- coding: utf-8 -*-
# __author__ = 'hao.gao'
# last_modify: 2019.09.24

import os,shutil,re
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from mbio.packages.bac_comp_genome.common_function import link_dir,link_file

class BacGbkModule(Module):
    """
    多个样品的gbk生成
    """

    def __init__(self, work_id):
        super(BacGbkModule, self).__init__(work_id)
        option = [
            {"name": "gene_path", "type": "string"},  # string类型的基因预测dir路径
            {"name": "genome_path", "type": "string"},  # string类型的组装序列dir路径
            {"name": "sample_list", "type": "infile", "format": "sequence.profile_table"},#文件中有样品信息
            {"name": "gff_path", "type": "string"}, #上传的gff，没有处理的gff路径
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
        if not self.option("gene_path"):
            raise OptionError("必须提供gene_path的路径！")
        if not self.option("genome_path"):
            raise OptionError("必须提供genome_path的路径！")
        return True

    def run_all_gbk(self):
        self.samples = self.get_sample(self.option("sample_list").prop['path'])
        for sample in sorted(self.samples.keys()):
            if self.samples[sample] in ['gff','Gff']:
                gbk = self.add_tool('bac_comp_genome.bac_gbk')
                opts = {
                    "gene_gff": self.option("gene_path") + sample + "/" + sample + "_CDS.gff",
                    "rrna_gff": self.option("gene_path") + sample + "/" + sample + "_rRNA.gff",
                    "trna_gff": self.option("gene_path") + sample + "/" + sample + "_tRNA.gff",
                    "genome_fa": self.option("genome_path") + sample + ".fna",
                    "pro_fa": self.option("gene_path") + sample + "/" + sample + "_CDS.faa",
                    "sample_name": sample,
                    "path": self.option("gff_path") + "/" + sample,
                    "type": "genome",
                }
                gbk.set_options(opts)
                self.all_modules.append(gbk)
            else:
                gbk = self.add_tool('bac_comp_genome.bac_gbk')
                opts = {
                    "gene_gff": self.option("gene_path") + sample + "/" + sample + "_CDS.gff",
                    "rrna_gff": self.option("gene_path") + sample + "/" + sample + "_rRNA.gff",
                    "trna_gff": self.option("gene_path") + sample + "/" + sample + "_tRNA.gff",
                    "genome_fa": self.option("genome_path") + sample + ".fna",
                    "pro_fa": self.option("gene_path") + sample + "/" + sample + "_CDS.faa",
                    "sample_name": sample,
                    "type": "bac",
                }
                gbk.set_options(opts)
                self.all_modules.append(gbk)
        if len(self.all_modules) >=1:
            self.on_rely(self.all_modules, self.set_output)
        else:
            self.all_modules[0].on("end", self.set_output)
        for module in self.all_modules:
            module.run()

    def set_output(self):
        for module in self.all_modules:
            if len(os.listdir(module.output_dir)) >1:
                link_dir(module.output_dir, self.output_dir)
        self.end()

    def get_sample(self,file):
        list = {}
        with open(file, "r") as f:
            lines = f.readlines()
            for line in lines[1:]:
                lin = line.strip().split("\t")
                list[lin[0]] = lin[1]
        return list

    def run(self):
        super(BacGbkModule, self).run()
        self.run_all_gbk()

    def end(self):
        super(BacGbkModule, self).end()
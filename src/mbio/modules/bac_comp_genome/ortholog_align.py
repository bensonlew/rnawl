# -*- coding: utf-8 -*-
# __author__ = 'hao.gao'
# last_modify: 2019.10.19

import os,shutil
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from mbio.packages.bac_comp_genome.common_function import link_dir,link_file

class OrthologAlignModule(Module):
    """
    根据不同方法构建基因树
    """

    def __init__(self, work_id):
        super(OrthologAlignModule, self).__init__(work_id)
        option = [
            {"name": "type", "type": "string", "default": "port"},  # 核酸nul or 蛋白port
            {"name": "homolog", "type": "infile", "format": "sequence.profile_table"}, #同源蛋白聚类cluster表
            {"name": "gene_dir", "type": "string"},
            {"name": "sample_list", "type": "string"},  # 文件中有样品信息
            {"name": "out", "type": "outfile", "format": "sequence.fasta"},  # 比对对齐的序列文件
        ]
        self.add_option(option)
        self.get_fasta = self.add_tool("bac_comp_genome.get_ortholog")
        self.merge_fasta = self.add_tool("bac_comp_genome.merge_fasta")
        self.all_modules = []

    def check_options(self):
        """
        检查参数
        :return:
        """
        if not self.option("sample_list"):
            raise OptionError("必须提供样品的信息文件！")
        if not self.option("homolog").is_set:
            raise OptionError("必须提供homolog的文件！")
        return True

    def run_get_fasta(self):
        self.get_fasta.set_options({
            "type": self.option("type"),
            "homolog": self.option("homolog"),
            "gene_path": self.option("gene_dir"),
            "sample_list": self.option("sample_list"),
        })
        self.get_fasta.on("end", self.run_align)
        self.get_fasta.run()

    def run_align(self):
        files = os.listdir(self.get_fasta.option("out").prop['path'])
        for file in files:
            self.align = self.add_tool("bac_comp_genome.align_triml")
            opts = {
                "seq": self.get_fasta.option("out").prop['path'] + "/" + file,
            }
            self.align.set_options(opts)
            self.all_modules.append(self.align)
        if len(self.all_modules) >1:
            self.on_rely(self.all_modules, self.run_merge_fasta)
        else:
            self.all_modules[0].on("end",self.set_output)
        for module in self.all_modules:
            module.run()

    def run_merge_fasta(self):
        if os.path.exists(self.work_dir + "/merge"):
           shutil.rmtree(self.work_dir + "/merge")
        os.mkdir(self.work_dir + "/merge")
        n =0
        for module in self.all_modules:
            n +=1
            link_file(module.option("out").prop["path"], self.work_dir + "/merge/align" + str(n)+".fasta")
        self.merge_fasta.set_options({
            "fasta_dir": self.work_dir + "/merge",
            "sample_list": self.option("sample_list"),
        })
        self.merge_fasta.on("end", self.set_output)
        self.merge_fasta.run()

    def set_output(self):
        if len(self.all_modules) > 1:
            link_file(self.merge_fasta.output_dir + "/all.align.fasta", self.output_dir + "/all.align.fasta")
            self.option("out", self.output_dir + "/all.align.fasta")
        else:
            link_file(self.all_modules[0].output_dir + "/all.align_last.fna", self.output_dir + "/all.align.fasta")
            self.option("out", self.output_dir + "/all.align.fasta")
        self.end()

    def run(self):
        super(OrthologAlignModule, self).run()
        self.run_get_fasta()

    def end(self):
        super(OrthologAlignModule, self).end()
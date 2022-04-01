# -*- coding: utf-8 -*-
# __author__ = 'hao.gao'
# last_modify: 2019.10.15

import os,shutil
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from mbio.packages.bac_comp_genome.common_function import link_dir,link_file

class SpeciesTreeModule(Module):
    """
    根据不同方法构建物种树
    """

    def __init__(self, work_id):
        super(SpeciesTreeModule, self).__init__(work_id)
        option = [
            {"name": "analysis_type", "type": "string", "default": "s16"},  # s16、ortholog、hous_keeping
            {"name": "type", "type": "string", "default": "port"},  # 核酸nul or 蛋白port
            {"name": "homolog", "type": "infile", "format": "sequence.profile_table"}, #同源蛋白聚类cluster表
            {"name": "gene_path", "type": "string"},
            {"name": "sample_list", "type": "string"},  # 文件中有样品信息
            {"name": "16s_dir", "type": "infile", "format": "sequence.fasta_dir"},  # 16s的文件夹
            {"name": "core_gene", "type": "infile", "format": "sequence.fasta"},  # core_gene的合并序列
            {"name": "method", "type": "string", "default": "NJ"},  # NJ or ML
            {"name": "out_group", "type": "string"},  # 外类群在进化文件中名称
            {"name": "bootstrap", "type": "int", "default": 1000},  # 外类群在进化文件中名称
            {"name": "merit", "type": "string", "default": "BIC"},  # 模型选择评估标准，BIC、AIC和AICc
        ]
        self.add_option(option)

    def check_options(self):
        """
        检查参数
        :return:
        """
        pass

    def run_s16(self):
        self.s16 = self.add_module("bac_comp_genome.rrna_tree")
        opts = {
            "16s_dir": self.option("16s_dir"),
            "method": self.option("method"),
            "bootstrap": self.option("bootstrap"),
        }
        if self.option("method") in ["ML"]:
            opts["merit"] = self.option("merit")
            if self.option("out_group"):
                opts["out_group"] = self.option("out_group")
        self.s16.set_options(opts)
        self.s16.on("end", self.set_output)
        self.s16.run()

    def run_get_ortholog(self):
        self.get_ortholog_fasta = self.add_module("bac_comp_genome.ortholog_align")
        opts = {
            "type": self.option("type"),
            "homolog": self.option("homolog"),
            "gene_dir": self.option("gene_path"),
            "sample_list": self.option("sample_list"),
        }
        self.get_ortholog_fasta.set_options(opts)
        self.get_ortholog_fasta.on("end", self.run_tree)
        self.get_ortholog_fasta.run()

    def run_align(self):
        self.align = self.add_tool("bac_comp_genome.align_triml")
        seq = ''
        if self.option("analysis_type") in ["house_keeping"]:
            seq = self.option("core_gene")
        opts = {
            "seq": seq,
        }
        self.align.set_options(opts)
        self.align.on("end", self.run_tree)
        self.align.run()

    def run_tree(self):
        self.tree = self.add_tool("bac_comp_genome.iqtree")
        seq = ''
        if self.option("analysis_type") in ["house_keeping"]:
            seq = self.align.option("out")
        elif self.option("analysis_type") in ["ortholog"]:
            seq = self.get_ortholog_fasta.option("out")
        opts = {
            "fa": seq,
            "bootstrap": self.option("bootstrap"),
            "merit": self.option("merit"),
        }
        if self.option("out_group"):
            opts["out_group"] = self.option("out_group")
        self.tree.set_options(opts)
        self.tree.on("end", self.set_output)
        self.tree.run()

    def set_output(self):
        if self.option("analysis_type") in ["s16"]:
            link_dir(self.s16.output_dir, self.output_dir)
        elif self.option("analysis_type") in ["ortholog"]:
            link_dir(self.tree.output_dir, self.output_dir)
        elif self.option("analysis_type") in ["house_keeping"]:
            link_dir(self.tree.output_dir, self.output_dir)
        self.end()

    def run(self):
        super(SpeciesTreeModule, self).run()
        if self.option("analysis_type") in ["s16"]:
            self.run_s16()
        elif self.option("analysis_type") in ["ortholog"]:
            self.run_get_ortholog()
        elif self.option("analysis_type") in ["house_keeping"]:
            self.run_align()

    def end(self):
        super(SpeciesTreeModule, self).end()
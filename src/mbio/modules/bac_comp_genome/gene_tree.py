# -*- coding: utf-8 -*-
# __author__ = 'hao.gao'
# last_modify: 2019.10.15

import os,shutil
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from mbio.packages.bac_comp_genome.common_function import link_dir,link_file

class GeneTreeModule(Module):
    """
    根据不同方法构建基因树
    """

    def __init__(self, work_id):
        super(GeneTreeModule, self).__init__(work_id)
        option = [
            {"name": "analysis_type", "type": "string", "default": "homolog"},  # ortholog、house_keeping
            {"name": "type", "type": "string", "default": "port"},  # 核酸nul or 蛋白port
            {"name": "homolog", "type": "infile", "format": "sequence.profile_table"}, #同源蛋白聚类cluster表
            {"name": "core_gene_blast", "type": "infile", "format": "sequence.profile_table"},  # core_gene的合并blast结果表
            {"name": "gene_path", "type": "string"},
            {"name": "sample_list", "type": "string"},  # 文件中有样品信息
            {"name": "clusters", "type": "string"},  # 逗号隔开的cluster;如：Cluster1，cluster2
            {"name": "core_names", "type": "string"},  # 逗号隔开的core names;如：rpsB,rplE
            {"name": "method", "type": "string", "default": "ML"},  # NJ or ML
            {"name": "out_group", "type": "string"},  # 外类群在进化文件中名称
            {"name": "bootstrap", "type": "int", "default": 1000},  # 外类群在进化文件中名称
            {"name": "merit", "type": "string", "default": "BIC"},  # 模型选择评估标准，BIC、AIC和AICc
        ]
        self.add_option(option)
        self.all_modules = []

    def check_options(self):
        """
        检查参数
        :return:
        """
        if not self.option("sample_list"):
            raise OptionError("必须提供样品的信息文件！")
        if not self.option("analysis_type"):
            raise OptionError("必须提供analysis_type的方法！")
        return True

    def run_get_homolog(self):
        self.get_homolog_fasta = self.add_tool("bac_comp_genome.get_gene_homolog")
        opts = {
            "type": self.option("type"),
            "homolog": self.option("homolog"),
            "clusters": self.option("clusters"),
            "gene_path": self.option("gene_path"),
            "sample_list": self.option("sample_list"),
        }
        self.get_homolog_fasta.set_options(opts)
        self.get_homolog_fasta.on("end", self.run_align)
        self.get_homolog_fasta.run()

    def run_get_housekeeping(self):
        self.get_housekeeping_fasta = self.add_tool("bac_comp_genome.get_housekeeping")
        opts = {
            "type": self.option("type"),
            "core_gene_blast": self.option("core_gene_blast"),
            "core_names": self.option("core_names"),
            "sample_list": self.option("sample_list"),
        }
        if self.option("out_group"):
            opts['out_group'] = self.option("out_group")
        self.get_housekeeping_fasta.set_options(opts)
        self.get_housekeeping_fasta.on("end", self.run_align)
        self.get_housekeeping_fasta.run()

    def run_align(self):
        self.align = self.add_tool("bac_comp_genome.align_triml")
        seq = ''
        if self.option("analysis_type") in ["ortholog"]:
            seq = self.get_homolog_fasta.option("out")
        elif self.option("analysis_type") in ["house_keeping"]:
            seq = self.get_housekeeping_fasta.option("out")
        opts = {
            "seq": seq,
        }
        self.align.set_options(opts)
        self.align.on("end", self.run_tree)
        self.align.run()

    def run_tree(self):
        self.tree = self.add_tool("bac_comp_genome.iqtree")
        opts = {
            "fa": self.align.option("out"),
            "bootstrap": self.option("bootstrap"),
            "merit": self.option("merit"),
        }
        self.tree.set_options(opts)
        self.tree.on("end", self.set_output)
        self.tree.run()

    def set_output(self):
        if self.option("analysis_type") in ["ortholog"]:
            link_dir(self.tree.output_dir, self.output_dir)
        elif self.option("analysis_type") in ["house_keeping"]:
            link_dir(self.tree.output_dir, self.output_dir)
        self.end()

    def run(self):
        super(GeneTreeModule, self).run()
        if self.option("analysis_type") in ["ortholog"]:
            self.run_get_homolog()
        elif self.option("analysis_type") in ["house_keeping"]:
            self.run_get_housekeeping()

    def end(self):
        super(GeneTreeModule, self).end()
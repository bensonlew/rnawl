# -*- coding: utf-8 -*-
# __author__ = 'hao.gao'
# last_modify: 2019.10.14

import os,shutil
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from mbio.packages.bac_comp_genome.common_function import merge_16s,link_file,link_dir

class RrnaTreeModule(Module):
    """
    16s构建进化树
    """

    def __init__(self, work_id):
        super(RrnaTreeModule, self).__init__(work_id)
        option = [
            {"name": "16s_dir", "type": "infile", "format": "sequence.fasta_dir"},  # 16s的文件夹
            {"name": "method", "type": "string","default": "NJ"}, # NJ or ML
            {"name": "out_group", "type": "string"},  # 外类群在进化文件中名称
            {"name": "bootstrap", "type": "int", "default": 1000},  # 外类群在进化文件中名称
            {"name": "merit", "type": "string"},  # 模型选择评估标准，BIC、AIC和AICc
        ]
        self.add_option(option)
        self.align =self.add_tool("bac_comp_genome.ssu_align")
        self.all_modules = []

    def check_options(self):
        """
        检查参数
        :return:
        """
        if not self.option("16s_dir").is_set:
            raise OptionError("必须提供样品的16s_dir信息文件！")
        if not self.option("method"):
            raise OptionError("必须提供method的路径！")
        else:
            if self.option("method") not in ["NJ", "ML"]:
                raise OptionError("必须提供正确method的方法 NJ or ML！")
        return True

    def run_align(self):
        opts = {
            "16s_fa": self.work_dir + "/all.16s.fasta",
            }
        self.align.set_options(opts)
        self.align.run()

    def run_nj_tree(self):
        self.s16_nj = self.add_tool("graph.phy_tree")
        self.s16_nj.set_options({
            "align": "false",
            "fasta": self.align.option("out"),
        })
        self.s16_nj.on("end", self.set_output)
        self.s16_nj.run()

    def run_ml_tree(self):
        self.s16_ml = self.add_tool("bac_comp_genome.iqtree")
        opts = {
            "fa": self.align.option("out"),
            "bootstrap": self.option("bootstrap"),
            "merit": self.option("merit"),
             }
        if self.option("out_group"):
            opts["out_group"] = self.option("out_group")
        self.s16_ml.set_options(opts)
        self.s16_ml.on("end", self.set_output)
        self.s16_ml.run()

    def set_output(self):
        if self.option("method") in ["NJ"]:
            link_file(self.s16_nj.output_dir + "/phylo_tree.nwk", self.output_dir + "/all.tree.nwk")
        elif self.option("method") in ["ML"]:
            link_file(self.s16_ml.output_dir + "/all.tree.nwk", self.output_dir + "/all.tree.nwk")
            link_file(self.s16_ml.output_dir + "/all.assess_result.xls", self.output_dir + "/all.assess_result.xls")
        self.end()

    def run(self):
        super(RrnaTreeModule, self).run()
        num = merge_16s(self.option("16s_dir").prop['path'], self.work_dir + "/all.16s.fasta")
        if num >=4:
            if self.option("method") in ["NJ"]:
                self.align.on("end", self.run_nj_tree)
            elif self.option("method") in ["ML"]:
                self.align.on("end", self.run_ml_tree)
            self.run_align()
        else:
            self.end()

    def end(self):
        super(RrnaTreeModule, self).end()
# -*- coding: utf-8 -*-
# __author__ = 'hao.gao'
# last_modify: 2019.09.10

import os,shutil
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from mbio.packages.bac_comp_genome.common_function import link_file

class AnnoKeggModule(Module):
    """
    单个组细菌基因基kegg数据库注释
    """
    def __init__(self, work_id):
        super(AnnoKeggModule, self).__init__(work_id)
        option = [
            {"name": "query", "type": "infile", "format": "sequence.fasta"},  # 序列
            {"name": "sample", "type": "string", "default": ""},  # 样品名
            {"name": "kegg", "type": "outfile", "format": "sequence.profile_table"},
            {"name": "secretory", "type": "outfile", "format": "sequence.profile_table"},
        ]
        self.add_option(option)
        self.top_diamond = self.add_module('align.meta_diamond')
        self.kegg_anno = self.add_tool('bac_comp_genome.anno_kegg')
        self.secretory = self.add_tool("bac_comp_genome.anno_secretory")

    def check_options(self):
        """
        检查参数
        :return:
        """
        if not self.option("query").is_set:
            raise OptionError("必须提供基因序列文件")
        if self.option("sample") == "":
            raise OptionError("必须提供样品名称")
        return True

    def run_top_dimond(self):
        opts = {
            'query': self.option('query'),
            'query_type': 'prot',
            'database': "kegg",
            'outfmt': 5
        }
        self.top_diamond.set_options(opts)
        self.top_diamond.on("end", self.run_kegg_anno)
        self.top_diamond.run()

    def run_kegg_anno(self):
        opts = {
            'kegg_xml': self.top_diamond.option('outxml')
        }
        self.kegg_anno.set_options(opts)
        self.kegg_anno.on("end", self.run_secretory)
        self.kegg_anno.run()

    def run_secretory(self):
        opts = {
            "kegg_anno":self.kegg_anno.output_dir + "/gene_kegg_anno.xls"
        }
        self.secretory.set_options(opts)
        self.secretory.on("end", self.set_output)
        self.secretory.run()

    def run(self):
        super(AnnoKeggModule, self).run()
        self.run_top_dimond()

    def set_output(self, event):
        link_file(self.kegg_anno.output_dir + "/gene_kegg_anno.xls", self.output_dir + "/" +self.option('sample') + ".kegg_anno.xls")
        link_file(self.secretory.output_dir + "/secretory_system.xls", self.output_dir + "/" +self.option('sample') + ".secretory_system.xls")
        self.option("secretory", self.output_dir + "/" +self.option('sample') + ".secretory_system.xls")
        self.option("kegg", self.output_dir + "/" + self.option('sample') + ".kegg_anno.xls")
        self.end()

    def end(self):
        super(AnnoKeggModule, self).end()
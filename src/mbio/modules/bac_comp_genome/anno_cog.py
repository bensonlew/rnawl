# -*- coding: utf-8 -*-
# __author__ = 'hao.gao'
# last_modify: 2019.09.10

import os,shutil
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from mbio.packages.bac_comp_genome.common_function import link_file

class AnnoCogModule(Module):
    """
    单个细菌基因组的COG注释
    """
    def __init__(self, work_id):
        super(AnnoCogModule, self).__init__(work_id)
        option = [
            {"name": "query", "type": "infile", "format": "sequence.fasta"},  # 序列
            {"name": "sample", "type": "string", "default": ""},  # 样品名
            {"name": "cog", "type": "outfile", "format": "sequence.profile_table"},
        ]
        self.add_option(option)
        self.top_diamond = self.add_module('align.meta_diamond')
        self.cog_anno = self.add_tool('bac_comp_genome.anno_cog')

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
            'database': "eggnog",
            'outfmt': 5
        }
        self.top_diamond.set_options(opts)
        self.top_diamond.on("end", self.run_cog_anno)
        self.top_diamond.run()

    def run_cog_anno(self):
        opts = {
            'cog_xml': self.top_diamond.option('outxml')
        }
        self.cog_anno.set_options(opts)
        self.cog_anno.on("end", self.set_output)
        self.cog_anno.run()

    def run(self):
        super(AnnoCogModule, self).run()
        self.run_top_dimond()

    def set_output(self, event):
        link_file(self.cog_anno.output_dir + "/gene_cog_anno.xls", self.output_dir + "/" + self.option("sample") + ".cog_anno.xls")
        self.option("cog", self.output_dir + "/" + self.option("sample") + ".cog_anno.xls")
        self.end()

    def end(self):
        super(AnnoCogModule, self).end()
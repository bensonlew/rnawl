# -*- coding: utf-8 -*-
# __author__ = 'gao.hao'
# last_modify: 2019.09.23

import os,shutil
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from mbio.packages.bac_comp_genome.common_function import link_file

class AnnoCazyModule(Module):
    """
    单个细菌cazy注释
    """

    def __init__(self, work_id):
        super(AnnoCazyModule, self).__init__(work_id)
        option = [
            {"name": "query", "type": "infile", "format": "sequence.fasta"},  # 序列
            {"name": "sample", "type": "string", "default": ""},  # 样品名
            {"name": "cazy", "type": "outfile", "format": "sequence.profile_table"}
        ]
        self.add_option(option)
        self.hmmscan = self.add_module('align.hmmscan')
        self.cazy_anno = self.add_tool('bac_comp_genome.anno_cazy')

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

    def run_hmmscan(self):
        opts = {
            'query': self.option('query'),
            'database': 'cazy',
        }
        self.hmmscan.set_options(opts)
        self.hmmscan.on("end", self.run_cazy_anno)
        self.hmmscan.run()

    def run_cazy_anno(self):
        opts = {
            'hmmscan_result': self.hmmscan.option('align_result'),
            'add_score': 'True'
        }
        self.cazy_anno.set_options(opts)
        self.cazy_anno.on("end", self.set_output)
        self.cazy_anno.run()

    def run(self):
        super(AnnoCazyModule, self).run()
        self.run_hmmscan()

    def set_output(self):
        link_file(self.cazy_anno.output_dir + "/" + "gene_cazy_parse_anno.xls", self.output_dir + "/" + self.option("sample") + ".cazy_anno.xls")
        self.option("cazy", self.output_dir + "/" + self.option("sample") + ".cazy_anno.xls")
        self.end()

    def end(self):
        super(AnnoCazyModule, self).end()
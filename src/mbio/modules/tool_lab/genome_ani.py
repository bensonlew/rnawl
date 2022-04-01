# -*- coding: utf-8 -*-
# __author__ = 'hao.gao'
# last_modify: 2019.10.09

import os,shutil
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from mbio.packages.bac_comp_genome.common_function import link_dir,link_file,add_merge,add_fasta,chang_value
import pandas as pd

class GenomeAniModule(Module):
    """
    多个样品的core gene的获取
    """

    def __init__(self, work_id):
        super(GenomeAniModule, self).__init__(work_id)
        option = [
            {"name": "seq_dir", "type": "infile", "format": "sequence.fasta_dir"},  # 输入参考序列文件夹
            {"name": "out", "type": "outfile", "format": "sequence.profile_table"},  #输出ani计算结果
            {"name": "method", "type": "string", "default": "ANIm"}, #ANIm,ANIb,ANIblastall,TETRA
            {"name": "linkage", "type": "string", "default": "average"},
        ]
        self.add_option(option)
        self.ani = self.add_tool("bac_comp_genome.genome_ani")
        self.all_modules = []

    def check_options(self):
        """
        检查参数
        :return:
        """
        if not self.option("seq_dir").is_set:
            raise OptionError("必须提供seq_dir文件夹！")
        return True

    def run_ani(self):
        self.ani.set_options({
            "seq_dir": self.option("seq_dir"),
            "method":self.option("method")
        })
        self.ani.on("end", self.run_hcluster)
        self.ani.run()

    def run_hcluster(self):
        self.hcluster = self.add_tool('bac_comp_genome.hcluster')
        if self.option("method") in ['TETRA']:
            chang_value(self.ani.option('out').prop['path'], self.work_dir + "/all.dismisily.xls", type=1)
        else:
            chang_value(self.ani.option('out').prop['path'], self.work_dir + "/all.dismisily.xls")
        self.logger.info(self.ani.option('out').prop['path'])
        self.hcluster.set_options({
            'dis_matrix': self.work_dir + "/all.dismisily.xls",
            'linkage': self.option("linkage")
        })
        self.hcluster.on('end', self.set_output)
        self.hcluster.run()

    def set_output(self):
        self.get_file(self.ani.output_dir + "/all.ANI_summary.xls", self.work_dir + "/all.ANI_summary.xls")
        link_file(self.work_dir + "/all.ANI_summary.xls", self.output_dir + "/all.ANI_summary.xls")
        link_file(self.hcluster.output_dir + "/hcluster.tre", self.output_dir + '/hcluster.tre')
        self.end()

    def get_file(self, file, out):
        table = pd.read_table(file, sep='\t', header=0)
        table = table.rename(columns={'Unnamed: 0': 'samples'})
        table.to_csv(out, sep='\t', header=True, index=False)

    def run(self):
        super(GenomeAniModule, self).run()
        self.run_ani()

    def end(self):
        super(GenomeAniModule, self).end()
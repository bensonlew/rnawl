# -*- coding: utf-8 -*-
# __author__ = 'hao.gao'
# last_modify: 2019.10.09

import os,shutil
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from mbio.packages.bac_comp_genome.common_function import link_dir,link_file,add_merge,add_fasta,chang_value

class GenomeAaiModule(Module):
    """
    多个样品的core gene的获取
    """

    def __init__(self, work_id):
        super(GenomeAaiModule, self).__init__(work_id)
        option = [
            {"name": "seq_dir", "type": "infile", "format": "sequence.fasta_dir"},  # 输入参考序列文件夹
            {"name": "evalue", "type": "string", "default": "1e-3"},
            {"name": "identity", "type": 'int', "default": 30},
            {"name": "aln_len", "type": 'int', "default": 70},
            {"name": "file_ext", "type": "string", "default": "fna"},
            {"name": "linkage", "type": "string", "default": "average"},
        ]
        self.add_option(option)
        self.aai = self.add_tool("bac_comp_genome.genome_aai")
        self.all_modules = []

    def check_options(self):
        """
        检查参数
        :return:
        """
        if not self.option("seq_dir").is_set:
            raise OptionError("必须提供seq_dir文件夹！")
        return True

    def run_aai(self):
        self.aai.set_options({
            "seq_dir": self.option("seq_dir"),
            "evalue": self.option("evalue"),
            "identity": self.option("identity"),
            "aln_len": self.option("aln_len"),
            "file_ext": self.option("file_ext"),
        })
        self.aai.on("end", self.run_hcluster)
        self.aai.run()

    def run_hcluster(self):
        self.hcluster = self.add_tool('bac_comp_genome.hcluster')
        chang_value(self.aai.option('out').prop['path'], self.work_dir + "/all.dismisily.xls")
        self.logger.info(self.aai.option('out').prop['path'])
        self.hcluster.set_options({
            'dis_matrix': self.work_dir + "/all.dismisily.xls",
            'linkage': self.option("linkage")
        })
        self.hcluster.on('end', self.set_output)
        self.hcluster.run()

    def set_output(self):
        link_file(self.aai.output_dir + "/aai_summary.xls", self.output_dir + "/aai_summary.xls")
        link_file(self.hcluster.output_dir + "/hcluster.tre", self.output_dir + '/hcluster.tre')
        self.end()

    def run(self):
        super(GenomeAaiModule, self).run()
        self.run_aai()

    def end(self):
        super(GenomeAaiModule, self).end()
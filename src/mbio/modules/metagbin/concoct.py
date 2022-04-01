# -*- coding: utf-8 -*-
# __author__ = 'gaohao'

import os
import re
import time, shutil
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from mbio.packages.metagbin.common_function import link_dir


class ConcoctModule(Module):
    """
    按样品合并bam文件
    author: gaohao
    last_modify: 2019.07.029
    """
    def __init__(self, work_id):
        super(ConcoctModule, self).__init__(work_id)
        options = [
            {"name": "bam_dir", "type": "infile", "format": "metagbin.bam_dir"},  # bam的文件夹
            {'name': 'specimen_info', 'type': 'infile', 'format': 'meta_genomic.specimen_info'},
            {"name": "minContig", "type": "string", "default": "1000"},  # metabat2最小contigs
            {"name": "contig_fa", "type": "infile", "format": "sequence.fasta"},  # metabat2输入文件contigs.fa
            {"name": "concoct_bin", "type": "outfile", "format": "sequence.fasta_dir"},  # 生成bin的目录
        ]
        self.modules = []
        self.add_option(options)
        self.bam = self.add_module('metagbin.sample_merge_bam')
        self.concoct = self.add_tool('metagbin.concoct')

    def check_options(self):
        """
        检查参数
        :return:
        """
        if not self.option('bam_dir').is_set:
            raise OptionError('必须输入bam_dir文件夹')

    def run_bam(self):
        self.bam.set_options({
            'bam_dir': self.option("bam_dir"),
            'specimen_info':self.option("specimen_info"),
        })
        self.bam.on("end",self.run_concoct)
        self.bam.run()

    def run_concoct(self):
        opts = {
            "minContig": self.option("minContig"),
            "contig_fa": self.option("contig_fa"),
            "bam_dir": self.bam.output_dir + "/bam_sort",
        }
        self.concoct.set_options(opts)
        self.concoct.on("end",self.set_output)
        self.concoct.run()

    def run(self):
        """
        运行
        :return:
        """
        super(ConcoctModule, self).run()
        self.run_bam()

    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        if os.path.exists(self.output_dir + '/concoct_bin'):
            shutil.rmtree(self.output_dir + '/concoct_bin')
        link_dir(self.concoct.output_dir + "/concoct_bin", self.output_dir + '/concoct_bin')
        self.option('concoct_bin', self.output_dir + '/concoct_bin')
        self.end()

    def end(self):
        super(ConcoctModule, self).end()
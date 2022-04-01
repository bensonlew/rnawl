#-*- coding: utf-8 -*-
# __author__ = 'gaohao'

import os
import re,shutil
import time
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from mbio.packages.metagbin.common_function import link_dir

class GetBinCoregeneModule(Module):
    """
    宏基因组binning获取bin的core_gene
    """
    def __init__(self, work_id):
        super(GetBinCoregeneModule, self).__init__(work_id)
        options = [
            {"name": "scaftig", "type": "infile", "format": "sequence.fasta"},
            {"name": "sample_name", "type": "string"},  # 样本的名称
            {"name": "metabat_depth", "type": "outfile", "format": "sequence.profile_table"},
        ]
        self.add_option(options)
        self.metagene = self.add_tool('metagbin.metagene')
        self.coregene = self.add_tool('metagbin.get_core_gene')

    def check_options(self):
        """
        重写参数检测函数
        """
        if not self.option('scaftig').is_set:
            raise OptionError('必须输入scaftig序列文件！')
        if not self.option('sample_name'):
            raise OptionError('必须输入sample_name参数！')

    def run_metagene(self):
        """
        进行基因预测
        :return:
        """
        opts = ({
            'cut_more_scaftig': self.option('scaftig'),
            'sample_name': self.option('sample_name')
        })
        self.metagene.set_options(opts)
        self.metagene.on('end', self.run_get_coregene)
        self.metagene.run()

    def run_get_coregene(self):
        """
        获取看家基因
        :return:
        """
        opts = ({
            'seq_faa': self.metagene.option('faa'),
            "seq_gff": self.metagene.option('gff'),
            'sample_name':self.option('sample_name'),
        })
        self.coregene.set_options(opts)
        self.coregene.on("end",self.set_output)
        self.coregene.run()

    def set_output(self):
        """
        设置结果目录
        :return:
        """
        link_dir(self.coregene.output_dir,self.output_dir)
        self.end()

    def run(self):
        super(GetBinCoregeneModule, self).run()
        self.run_metagene()

    def end(self):
        super(GetBinCoregeneModule, self).end()
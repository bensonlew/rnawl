#-*- coding: utf-8 -*-
# __author__ = 'gaohao'

import os
import re,shutil
import time
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from mbio.packages.metagbin.common_function import link_dir

class MetagbinBinStatModule(Module):
    """
    宏基因组binning 结果注释、统计
    """
    def __init__(self, work_id):
        super(MetagbinBinStatModule, self).__init__(work_id)
        options = [
            {"name": "bin_dir", "type": "infile", "format": "sequence.fasta_dir"},  # bin的文件目录
            {"name": "input_genome", "type": "infile", "format": "sequence.fasta"},  ##组装序列
            {"name": "metabat_depth", "type": "outfile", "format": "sequence.profile_table"},
        ]
        self.add_option(options)
        self.split_fastq_tools = []
        self.amphora = self.add_module('metagbin.amphora')
        self.checkm = self.add_module('metagbin.checkm')
        self.rrna = self.add_tool('metagbin.rrna_predict')
        self.step.add_steps('checkm', 'amphora', 'rrna')

    def check_options(self):
        """
        重写参数检测函数
        """
        if not self.option('bin_dir').is_set:
            raise OptionError('必须输入bin_dir序列文件夹')
        if not self.option('input_genome').is_set:
            raise OptionError('必须输入input_genome序列文件')

    def run_checkm(self):
        """
        对所有bins进行评估
        :return:
        """
        opts = ({
            'bin_dir': self.option('bin_dir'),
        })
        self.checkm.set_options(opts)
        self.checkm.on('end', self.set_output, "checkm")
        self.checkm.run()

    def run_amphora(self):
        """
        进行物种注释
        :return:
        """
        opts = ({
            'bin_dir': self.option('bin_dir'),
            'taxon_file':self.checkm.output_dir + "/all.bin.summary.xls",
        })
        self.amphora.set_options(opts)
        self.amphora.on('end', self.set_output, "amphora")
        self.amphora.run()

    def run_rrna(self):
        """
        进行rRNA预测，提取16s
        :return:
        """
        opts = ({
            'input_genome': self.option('input_genome'),
            'metabat_depth':self.option('metabat_depth'),
            'bin_dir': self.option('bin_dir'),
        })
        self.rrna.set_options(opts)
        self.rrna.on('end', self.set_output, "rrna")
        self.rrna.run()

    def set_output(self,event):
        """
        设置结果目录
        :return:
        """
        obj = event['bind_object']
        if event['data'] == 'checkm':
            for i in ['all.marker.xls','all.bin.summary.xls']:
                if os.path.exists(self.output_dir + '/' + i):
                    os.remove(self.output_dir + '/' + i)
                os.link(obj.output_dir + '/' + i,self.output_dir + '/' + i)
        if event['data'] == 'amphora':
            for i in ['summary.anno.xls']:
                if os.path.exists(self.output_dir + '/' + i):
                    os.remove(self.output_dir + '/' + i)
                os.link(obj.output_dir + '/' + i,self.output_dir + '/' + i)
        if event['data'] == 'rrna':
            for i in ['bins.rRNA.gff','bins.rRNA.fnn','all_bins.scf.xls','all_bins.16s.xls']:
                if os.path.exists(self.output_dir + '/' + i):
                    os.remove(self.output_dir + '/' + i)
                os.link(obj.output_dir + '/' + i,self.output_dir + '/' + i)
            self.end()

    def run(self):
        super(MetagbinBinStatModule, self).run()
        self.checkm.on('end',self.run_amphora)
        self.amphora.on('end', self.run_rrna)
        self.run_checkm()

    def end(self):
        super(MetagbinBinStatModule, self).end()
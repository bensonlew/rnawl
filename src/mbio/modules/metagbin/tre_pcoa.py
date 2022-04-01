#-*- coding: utf-8 -*-
# __author__ = 'gaohao'

import os
import re,shutil
import time
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from mbio.packages.metagbin.common_function import link_dir

class TrePcoaModule(Module):
    """
    宏基因组四核苷酸频率的pcoa分析
    """
    def __init__(self, work_id):
        super(TrePcoaModule, self).__init__(work_id)
        options = [
            {"name": "input_genome", "type": "infile", "format": "sequence.fasta"},
            {"name": "dis_method", "type": "string", "default": "bray_curtis"},
            {"name": "ellipse", "type": "string", "default": "F"},
        ]
        self.add_option(options)
        self.matrix = self.add_tool('meta.beta_diversity.distance_calc')
        self.checkm = self.add_tool('metagbin.checkm_tra')
        self.pcoa = self.add_tool('meta.beta_diversity.pcoa')
        self.step.add_steps("checkm","matrix","pcoa")

    def check_options(self):
        """
        重写参数检测函数
        """
        if not self.option('input_genome').is_set:
            raise OptionError('必须输入input_genome序列文件!')

    def run_checkm(self):
        """
        计算四核苷酸频率表，并倒置
        :return:
        """
        opts = ({
            'input_genome': self.option('input_genome'),
        })
        self.checkm.set_options(opts)
        self.checkm.on('end', self.set_output, "checkm")
        self.checkm.run()

    def matrix_run(self):
        """
        运行计算距离矩阵
        :return:
        """
        self.matrix.set_options({'method': self.option('dis_method'),
                                     'otutable': self.checkm.option("out")})
        self.matrix.on('end', self.set_output, 'matrix')
        self.matrix.run()

    def run_pcoa(self):
        """
        计算pcoa
        :return:
        """
        pcoa_options = {'dis_matrix': self.matrix.option('dis_matrix'), 'ellipse': self.option('ellipse')}
        self.pcoa.set_options(pcoa_options)
        self.pcoa.on('end', self.set_output, 'pcoa')
        self.pcoa.run()

    def set_output(self,event):
        """
        设置结果目录
        :return:
        """
        obj = event['bind_object']
        if event['data'] == 'checkm':
            for i in ['tetra.summary.xls']:
                if os.path.exists(self.output_dir + '/' + i):
                    os.remove(self.output_dir + '/' + i)
                os.link(obj.output_dir + '/' + i,self.output_dir + '/' + i)
        if event['data'] == 'matrix':
            for i in ['bray_curtis_tetra.pcoa.xls']:
                if os.path.exists(self.output_dir + '/' + i):
                    os.remove(self.output_dir + '/' + i)
                os.link(obj.output_dir + '/' + i,self.output_dir + '/' + i)
        if event['data'] == 'pcoa':
            for i in ['pcoa_sites.xls','pcoa_eigenvaluespre.xls', 'pcoa_eigenvalues.xls']:
                if os.path.exists(self.output_dir + '/' + i):
                    os.remove(self.output_dir + '/' + i)
                os.link(obj.output_dir + '/' + i,self.output_dir + '/' + i)

    def run(self):
        super(TrePcoaModule, self).run()
        self.pcoa.on('end',self.end)
        self.matrix.on('end', self.run_pcoa)
        self.checkm.on('end', self.matrix_run)
        self.run_checkm()

    def end(self):
        super(TrePcoaModule, self).end()
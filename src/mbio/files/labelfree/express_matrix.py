# -*- coding: utf-8 -*-
import os
import pandas as pd
from biocluster.iofile import File
from biocluster.core.exceptions import FileError
__author__ = "gdq"


class ExpressMatrixFile(File):
    def __init__(self):
        super(ExpressMatrixFile, self).__init__()

    def get_info(self):
        super(ExpressMatrixFile, self).get_info()
        gene_number, sample_number = self.parse_file()
        self.set_property('sample_number', sample_number)
        self.set_property('gene_number', gene_number)

    def parse_file(self):
        exp_file = self.prop['path']
        if not os.path.exists(exp_file):
            raise FileError('文件不存在：{}'.format(exp_file))
        table = pd.read_table(exp_file, index_col=0, header=0)
        gene_number = table.shape[0]
        if gene_number <= 2:
            raise FileError('表达矩阵近乎为空：{}'.format(exp_file))
        sample_number = table.shape[1]
        if sample_number <= 1:
            raise FileError('表达矩阵包含的样本太少，无法做差异分析: {}'.format(exp_file))
        return gene_number, sample_number

    def check(self):
        super(ExpressMatrixFile, self).check()
        self.get_info()

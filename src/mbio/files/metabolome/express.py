# -*- coding: utf-8 -*-
import os
import pandas as pd
from biocluster.iofile import File
from biocluster.core.exceptions import FileError


class ExpressFile(File):
    def __init__(self):
        super(ExpressFile, self).__init__()

    def get_info(self):
        super(ExpressFile, self).get_info()
        gene_number, sample_number = self.parse_file()
        self.set_property('sample_number', sample_number)
        self.set_property('gene_number', gene_number)

    def parse_file(self):
        exp_file = self.prop['path']
        if not os.path.exists(exp_file):
            raise FileError('文件不存在：%s', variables=(exp_file), code="44700501")
        table = pd.read_csv(exp_file, index_col=0, header=0, sep='\t')
        gene_number = table.shape[0]
        if gene_number < 2:
            raise FileError('表达矩阵只有一行：%s', variables=(exp_file), code="44700502")
        sample_number = table.shape[1]
        return gene_number, sample_number

    def check(self):
        super(ExpressFile, self).check()
        self.get_info()

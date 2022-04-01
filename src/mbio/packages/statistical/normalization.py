# -*- coding: utf-8 -*-
# __author__ = 'xieshichang'
import pandas as pd
import numpy as np


class Normalization(object):
    def __init__(self, table_in, by='row', norm_meth='relative', out=None):
        """二维数据表标准化
        Args:
            table_in (file): 输入二维数据表文件名，包含行列名
            by (str, optional): 标准化方向 row, col, both. Defaults to 'row'.
            norm_meth (str, optional): 标准化方法 relative, min-max, log10, z. Defaults to 'relative'.
            out (str, optional): 输出文件名称.
        """
        if isinstance(table_in, pd.DataFrame):
            self.table = table_in
        else:
            self.table = pd.read_csv(table_in, sep='\t', index_col=0)
        self.check_format()
        self.by = by
        self.table_out = out
        self.norm_meth = norm_meth
        self.method_mapping = {  # 不同标准化方法对应的处理函数，均是*按行*
            'relative': self.meth_rl,
            'minmax': self.meth_mm,
            'log': self.meth_log,
            'z': self.meth_z,
        }

    def run(self, outfile=False):
        if self.by in ['col', 'both']:
            self.table = self.table.T
        try:
            self.method_mapping[self.norm_meth]()
            if self.by == 'both':
                self.table = self.table.T
                self.method_mapping[self.norm_meth]()
        except Exception as e:
            raise Exception('标准化遇到错误{}'.format(e))
        if self.by == 'col':
            self.table = self.table.T
        self.table = self.table.dropna(how='any')
        if self.table_out:
            self.table.to_csv(self.table_out, sep='\t')
            return self.table_out
        else:
            return self.table

    def check_format(self):
        cols = self.table.columns
        errs = filter(lambda x: self.table[x].dtype not in [int, float], cols)
        if errs:
            raise Exception('输入文件的列[{}]中存在不是数字的项!'.format(errs))

    def meth_z(self):
        '''
        z-score 标准化方法 (x - mean) / std
        '''
        std = self.table.std(axis=1)
        mean = self.table.mean(axis=1)
        self.table = ((self.table.T - mean) / std).T.dropna(how='all')

    def meth_log(self):
        '''
        log10 标准化 log10(x) / log10(max)
        '''
        table_log10 = np.log10(self.table)
        t_max = self.table.max(axis=1)
        max_log10 = np.log10(t_max)
        self.table = (table_log10.T / max_log10).T.dropna(how='all')
        self.table[abs(self.table) == np.inf] = 0

    def meth_rl(self):
        '''
        relative 标准化 x / sum
        '''
        table_sum = self.table.sum(axis=1)
        self.table = (self.table.T / table_sum).T.dropna(how='all')

    def meth_mm(self):
        '''
        min-max 标准化 (x - min) / (max - min)
        '''
        table_min = self.table.min(axis=1)
        table_max = self.table.max(axis=1)
        min_max = table_max - table_min
        self.table = ((self.table.T - table_min) / min_max).T.dropna(how='all')

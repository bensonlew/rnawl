# -*- coding: utf-8 -*-
"""
@time    : 2018/10/25 11:16
@file    : diff_anno_matrix.py
@author  : zhipeng.zhao
@contact: 757049042@qq.com
"""
import pandas as pd
from biocluster.iofile import File
# from src.biocluster.iofile import File

class DiffAnnoMatrixFile(File):
    """
    差异基因列表文件
    """
    def __init__(self):
        super(DiffAnnoMatrixFile, self).__init__()

    def get_info(self):
        return super(DiffAnnoMatrixFile, self).get_info()


    def check(self):
        if super(DiffAnnoMatrixFile, self).check():
            self.get_info()
        # 检测文件是否为空
        with self.get_reader() as infile:
            lines = []
            flag = 0
            for line in infile:
                flag += 1
                if flag == 5:
                    break
                elif line:
                    lines.append(line)
        if len(lines) <= 1:
            return False
        return True

    @property
    def anno_matrix(self):
        gene_list = []
        df = pd.read_table(self.get_reader(), sep='\t', header=0, index_col=0)
        df = df.dropna(how='all')
        return df



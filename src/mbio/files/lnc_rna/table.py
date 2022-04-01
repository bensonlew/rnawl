# -*- coding: utf-8 -*-
# __author__ = "liubinxu"
import os

from biocluster.iofile import File
from mbio.files.lnc_rna.common import CommonFile

class TableFile(CommonFile):
    def __init__(self):
        super(TableFile, self).__init__()

    def check(self):
        super(TableFile, self).check()


    def choose_by_list(self, line_num, choose_list, out, header=False):
        '''
        根据id list筛选表格文件
        '''
        with open(self.path, 'r') as f, open(out, 'w') as fo:
            if header:
                fo.write(f.readline())
            for line in f:
                if line.strip().split("\t")[line_num-1] in choose_list:
                    fo.write(line)

    def add_by_file(self, line_num, file_add, out, header=True):
        '''
        根据文件添加
        '''
        import pandas as pd
        exp_df = pd.read_table(self.path, header=0, index_col=0)
        exp_add = pd.read_table(file_add, header=None, index_col=0)
        exp_add2 = exp_add.rename(columns={1: 'seq_id', 2: 'rna_type', 3: 'is_new'})
        exp_add2.loc[:, 'is_new'] =  exp_add2.loc[:, 'is_new'].replace({'known': 'false', 'novel': 'true'})
        exp_pd3 = exp_add2[['rna_type', 'is_new']]
        exp_df =  exp_df.join(exp_pd3)
        exp_df.to_csv(out, sep="\t", header=True)

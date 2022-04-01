#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@time    : 2019/5/30 10:59
@file    : geneset_rnk.py
"""

import csv
import os
import re
import pandas as pd

from biocluster.iofile import File


class GenesetRnkFile(File):

    def check(self):
        super(GenesetRnkFile, self).check()
        for name, value in self.csv_reader():
            float(value.strip())
        df = pd.read_table(self.path, header=None, index_col=None)
        df.sort_values(1, ascending=False)
        out_path = os.path.join(os.path.dirname(self.path), 'sorted_' + os.path.basename(self.path))
        df.to_csv(out_path, sep='\t', header=False, index=False)
        self.set_path(out_path)
        return True

    def hard_link(self, new_path):
        if os.path.isdir(os.path.dirname(new_path)):
            cmd = 'ln {old} {new}'.format(old=self.path, new=new_path)
            return_code = os.system(cmd)
            if return_code not in (0, '0'):
                cmd = 'cp {old} {new}'.format(old=self.path, new=new_path)
                os.system(cmd)
        else:
            self.set_error('文件创建link失败，请确认文件路径是否正确', code="43703102")

    def basename(self):
        return os.path.basename(self.path)

    def csv_reader(self, remove_header=False):
        with self.get_reader() as in_handler:
            if remove_header:
                in_handler.readline()
            for item in csv.reader(in_handler, delimiter='\t'):
                yield item

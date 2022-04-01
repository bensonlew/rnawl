# -*- coding: utf-8 -*-

import csv
import json
import os
import pandas as pd

from biocluster.iofile import File

class LncCommonFile(File):
    def __init__(self):
        super(LncCommonFile, self).__init__()

    def check(self):
        super(LncCommonFile, self).check()

    def hard_link(self, new_path):
        if os.path.isdir(os.path.dirname(new_path)):
            cmd = 'ln {old} {new}'.format(old=self.path, new=new_path)
            return_code = os.system(cmd)
            if return_code not in (0, '0'):
                cmd = 'cp {old} {new}'.format(old=self.path, new=new_path)
                os.system(cmd)
        else:
            raise Exception('文件创建link失败，请确认文件路径是否正确')

    def basename(self):
        return os.path.basename(self.path)

    def csv_reader(self, remove_header=False):
        with self.get_reader() as in_handler:
            if remove_header:
                in_handler.readline()
            for item in csv.reader(in_handler, delimiter='\t'):
                yield item

    def csv_dict_reader(self):
        with self.get_reader() as in_handler:
            for dic in csv.DictReader(in_handler, delimiter='\t'):
                yield dic

    def json_reader(self):
        with self.get_reader() as in_handler:
            return json.load(in_handler)

    def dataframe(self, sep='\t', header=0, index_col=None):
        return pd.read_table(self.path, sep=sep, header=header, index_col=index_col)

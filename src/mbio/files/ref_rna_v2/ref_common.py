# -*- coding: utf-8 -*-

import csv
import json
import os
import pandas as pd

from biocluster.iofile import File

class RefCommonFile(File):

    def hard_link(self, new_path):
        if os.path.isdir(os.path.dirname(new_path)):
            cmd = 'ln {old} {new}'.format(old=self.path, new=new_path)
            return_code = os.system(cmd)
            if return_code not in (0, '0'):
                cmd = 'cp {old} {new}'.format(old=self.path, new=new_path)
                os.system(cmd)
        else:
            self.set_error('文件创建link失败，请确认文件路径是否正确', code="43703402")

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

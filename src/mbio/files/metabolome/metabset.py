# -*- coding: utf-8 -*-
# __author__ = 'xuting'
import re
from biocluster.iofile import File
from biocluster.core.exceptions import FileError
from collections import defaultdict
import pandas as pd


class MetabsetFile(File):
    """
    定义总览表结果文件
    检查表头名称，个数，顺序
    """
    def __init__(self):
        super(MetabsetFile, self).__init__()

    def get_info(self):
        """
        获取文件属性
        """
        super(MetabsetFile, self).get_info()
        info = self.get_file_info()
        self.set_property("set_len", len(info))
        self.set_property("set_list", info)
        # self.set_property("group_scheme", info[1])
        # self.set_property("is_empty", info[2])

    def get_file_info(self):
        """
        获取metab_table文件的信息
        """
        data = pd.read_csv(self.prop['path'], header=None, index_col=0, sep='\t')
        try:
            data_list = data.index.tolist()
        except:
            raise FileError("代谢集格式不正确，应为一列代谢物id", code="44700201")
        return data_list

    def format_check(self):
        pass

    def check(self):
        if super(MetabsetFile, self).check():
            self.get_info()

if __name__ == "__main__":
    g = MetabsetFile()
    g.set_path("example.group")
    g.get_info()
    g.sub_group("example.group.sub", ["g1", "g3", "g4"])

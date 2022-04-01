# -*- coding: utf-8 -*-
# __author__ = 'xuting'
import re
from biocluster.iofile import File
from biocluster.core.exceptions import FileError
from collections import defaultdict
import pandas as pd


class MulMetabsetFile(File):
    """
    定义总览表结果文件
    检查表头名称，个数，顺序
    """
    def __init__(self):
        super(MulMetabsetFile, self).__init__()
        self._intersection = True

    def get_info(self):
        """
        获取文件属性
        """
        super(MulMetabsetFile, self).get_info()
        info = self.get_file_info()
        self.set_property("is_one", info[0])
        self.set_property("is_intersection", info[1])
        self.set_property("set_name", info[2])
        # self.set_property("sample_number", len(info[0]))
        # self.set_property("sample", info[0])
        # self.set_property("group_scheme", info[1])
        # self.set_property("is_empty", info[2])

    def get_file_info(self):
        """
        获取metab_table文件的信息
        """
        set_name = []
        set_list = []
        only_one_set = True
        is_intersection = True
        with open(self.prop['path'], 'r') as f:
            lines = f.readlines()
            for line in lines:
                line = line.strip().split("\t")
                set_name.append(line[0])
                one_set_list = line[1].split(',')
                set_list.append(one_set_list)
        if len(set_list) == 1:
            only_one_set = True
        elif len(set_list) == 2:
            only_one_set = False
            inter = set(set_list[0]) & set(set_list[1])
            if len(inter) == 0:
                is_intersection = False
            else:
                is_intersection = True
        return only_one_set, is_intersection, set_name

    def format_check(self):
        pass

    def check(self):
        if super(MulMetabsetFile, self).check():
            self.get_info()

if __name__ == "__main__":
    g = MulMetabsetFile()
    g.set_path("example.group")
    g.get_info()
    g.sub_group("example.group.sub", ["g1", "g3", "g4"])

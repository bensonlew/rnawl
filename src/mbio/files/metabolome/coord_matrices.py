# -*- coding: utf-8 -*-
# __author__ = 'liulinmeng'
import re
from biocluster.iofile import File
from biocluster.core.exceptions import FileError
from collections import defaultdict
import pandas as pd


class CoordMatricesFile(File):
    """
    定义procrustes分析模块的输入文件格式
    检查表头名称，个数，顺序
    """
    def __init__(self):
        super(CoordMatricesFile, self).__init__()
        self.data = None

    def get_info(self):
        """
        获取文件属性
        """
        #super(MetabAbunFile, self).get_info()
        #info = self.get_file_info()
        #self.set_property("metab_list", info)
        #self.set_property("metab_num", len(info))
        # self.set_property("sample_number", len(info[0]))
        # self.set_property("sample", info[0])
        # self.set_property("group_scheme", info[1])
        # self.set_property("is_empty", info[2])

    def get_file_info(self):
        """
        获取文件的信息
        """
        #self.data = pd.read_table(self.prop['path'],index_col=["metab_id"])
        #if len(self.data) == 0:
        #    raise FileError("表格为空")
        #metab_list = self.data.index.tolist()
        #return metab_list

    def format_check(self):
        pass
        #if True in self.data.isnull().any().tolist():
         #   raise FileError("表格中有缺失值，需检查")

    def check(self):
        if super(CoordMatricesFile, self).check():
            self.get_info()
            self.format_check()


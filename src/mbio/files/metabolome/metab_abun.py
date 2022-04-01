# -*- coding: utf-8 -*-
# __author__ = 'xuting'
import re
from biocluster.iofile import File
from biocluster.core.exceptions import FileError
from collections import defaultdict
import pandas as pd


class MetabAbunFile(File):
    """
    定义上传的metab_table格式文件
    检查表头名称，个数，顺序
    """
    def __init__(self):
        super(MetabAbunFile, self).__init__()
        self.data = None

    def get_info(self):
        """
        获取文件属性
        """
        super(MetabAbunFile, self).get_info()
        info = self.get_file_info()
        self.set_property("metab_list", info)
        self.set_property("metab_num", len(info))
        # self.set_property("sample_number", len(info[0]))
        # self.set_property("sample", info[0])
        # self.set_property("group_scheme", info[1])
        # self.set_property("is_empty", info[2])

    def get_file_info(self):
        """
        获取metab_table文件的信息
        """
        self.data = pd.read_csv(self.prop['path'],index_col=["metab_id"], sep='\t')
        if len(self.data) == 0:
            raise FileError("表格为空", code="44700101")
        metab_list = self.data.index.tolist()
        return metab_list

    def format_check(self):
        if True in self.data.isnull().any().tolist():
            raise FileError("表格中有缺失值，需检查", code="44700102")

    def check(self):
        if super(MetabAbunFile, self).check():
            self.get_info()
            self.format_check()


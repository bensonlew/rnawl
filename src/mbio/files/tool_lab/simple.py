# -*- coding: utf-8 -*-
# __author__ = 'xuting'
import re
from biocluster.iofile import File
from biocluster.core.exceptions import FileError
from collections import defaultdict
import pandas as pd
import os


class SimpleFile(File):
    """
    定义总览表结果文件
    检查表头名称，个数，顺序
    """
    def __init__(self):
        super(SimpleFile, self).__init__()

    def get_info(self):
        """
        获取文件属性
        """
        super(SimpleFile, self).get_info()
        data = self.get_file_info()
        self.format_check(data)



    def get_file_info(self):
        """
        获取metab_table文件的信息
        """
        if os.path.getsize(self.prop['path']) == 0:
            raise  Exception('上传文件为空')
        data = pd.read_table(self.prop['path'], header=None, index_col=0)
        return data



    def format_check(self,data):
        if len(data) > 0:
            pass
        else:
            raise Exception('data is null')


    def check(self):
        if super(SimpleFile, self).check():
            self.get_info()

if __name__ == "__main__":
    g = SimpleFile()
    g.set_path("example.group")
    g.get_info()


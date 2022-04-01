# -*- coding: utf-8 -*-
# __author__ = 'zhaobinbin'
import os
from biocluster.core.exceptions import FileError
from biocluster.iofile import Directory


class CatlogDirFile(Directory):
    """
    定义文件夹格式
    """
    def __init__(self):
        super(CatlogDirFile, self).__init__()

    def check(self):
        """
        检测文件夹是否满足要求，不满足时触发FileError异常
        """
        if not os.path.isdir(self.prop['path']):
                raise FileError("不是一个文件夹", code="45500203")
        return False
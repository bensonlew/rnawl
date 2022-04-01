# -*- coding: utf-8 -*-
# __author__ = 'hongyu.chen'
"""xls格式文件"""
from biocluster.iofile import File


class XlsFile(File):
    """
    txt类
    """
    def __init__(self):
        super(XlsFile, self).__init__()

    def check(self):
        if super(XlsFile, self).check():
            return True
        else:
            raise FileError("文件格式错误")
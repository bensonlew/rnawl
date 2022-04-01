# -*- coding: utf-8 -*-
# __author__ = 'zhouxuan'

from biocluster.iofile import File
from biocluster.core.exceptions import FileError
import os


class ProfileFile(File):
    """
    定义宏基因的各种profile文件
    """
    def __init__(self):
        super(ProfileFile, self).__init__()

    def check_info(self):
        """
        获取判断文件夹内容
        """
        with open(self.prop['path'], 'r') as r:
            n = 0
            number = 0
            for line in r:
                n += 1
                if n == 1:
                    number = len(line.split("\t"))
                else:
                    if len(line.split("\t")) != number:
                        raise FileError('文件中信息不全，请检查', code="42701201")
        return True

    def check(self):
        """
        检测文件夹是否满足要求，不满足时触发FileError异常
        :return:
        """
        if super(ProfileFile, self).check():
            self.check_info()
            return True

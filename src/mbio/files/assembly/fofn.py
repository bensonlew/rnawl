# -*- coding: utf-8 -*-
# __author__ = 'gaohao'

from biocluster.iofile import File
from biocluster.core.exceptions import FileError
import subprocess
import os
import re


class FofnFile(File):
    """
    定义gff格式文件
    """

    def __init__(self):
        super(FofnFile, self).__init__()

    def check(self):
        """
        检测文件是否满足要求，发生错误时应该触发FileError异常
        :return:
        """
        if super(FofnFile, self).check():
            self.get_info()
            return True

    def get_info(self):
        super(FofnFile, self).get_info()
        if re.search('\.fofn',self.prop['path']):
            with open(self.prop["path"], "r") as f:
                lines =f.readlines()
                if len(lines) != 3:
                    raise FileError("h5文件格式错误，h5文件应有3个")
                line = f.readline()
                for line in f:
                    line = line.strip()
                    lst = line.split("\t")
                    # print lst
                    if len(lst) != 1:
                        raise FileError("文件格式错误，文件应只有1列")
                    if re.search('\.h5', lst[0]):
                        return True
                    else:
                        raise FileError("文件格式错误，文件应为.h5后缀的文件")
        else:
            raise FileError("文件格式错误，文件应为.fofn后缀的文件")
        return True



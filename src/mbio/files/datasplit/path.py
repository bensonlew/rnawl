# -*- coding: utf-8 -*-
# __author__ = 'wangzhaoyue'

"""针对非多样性barcode拆分格式文件类"""

from biocluster.iofile import File
from biocluster.core.exceptions import FileError
import re


class PathFile(File):
    """
   针对非多样性barcode拆分格式文件类
    """
    def __init__(self):
        super(PathFile, self).__init__()

    def check(self):
        """
        检测文件是否满足要求，发生错误时应该触发FileError异常
        :return:
        """
        if super(PathFile, self).check():
            with open(self.prop["path"], "r") as f:
                lines = f.readlines()
                for line in lines:
                    tmp_len = line.strip().split("\t")
                    # if len(tmp_len) != 1:
                    #     raise FileError("{}行格式不正确，因为每行仅一列，且为路径，请核实".format(line))
                    # if not re.search(r'/', line):
                    #     raise FileError("{}行格式不正确，应该是绝对路径，请核实".format(line))
        else:
            raise FileError("文件格式错误")

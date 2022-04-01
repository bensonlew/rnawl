# -*- coding: utf-8 -*-
# __author__ = 'zhouxuan'

from biocluster.iofile import File
from biocluster.core.exceptions import FileError
import re

class HmmscanTableFile(File):
    """
    定义宏基因的各种cazy比对处理后的结果文件
    """
    def __init__(self):
        super(HmmscanTableFile, self).__init__()

    def check_info(self):
        """
        获取判断文件夹内容
        """
        with open(self.prop['path'], 'r') as r:
            for line in r:
                if not line.startswith("#"):
                    line = line.strip().split(' ')
                    if len(line) != 0:
                        continue
                    else:
                        raise FileError('文件中信息不全，请检查', code="42800201")
        return True

    def check(self):
        """
        检测文件夹是否满足要求，不满足时触发FileError异常
        :return:
        """
        if super(HmmscanTableFile, self).check():
            self.check_info()
            return True

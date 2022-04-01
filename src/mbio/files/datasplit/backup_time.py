# -*- coding: utf-8 -*-
# __author__ = 'xuting'

import re
from biocluster.iofile import File
from biocluster.core.exceptions import FileError


class BackupTimeFile(File):
    """
    backup模块中备份拆分的数据所使用的年份和月份
    """
    def __init__(self):
        super(BackupTimeFile, self).__init__()

    def get_info(self):
        """
        获取文件信息
        """
        super(BackupTimeFile, self).get_info()
        with open(self.prop['path'], 'r') as r:
            for line in r:
                line = line.rstrip('\n')
                line = re.split('\t', line)
                if len(line) != 2:
                    raise FileError("文件格式错误")
                if line[0] != "year" and line[0] != 'month':
                    raise FileError("文件格式错误")
                if line[0] == "year":
                    self.set_property("year", line[1])
                if line[0] == "month":
                    self.set_property("month", line[1])

    def check(self):
        """
        检测文件是否满足要求,发生错误时应该触发FileError异常
        """
        if super(BackupTimeFile, self).check():
            self.get_info()
            return True

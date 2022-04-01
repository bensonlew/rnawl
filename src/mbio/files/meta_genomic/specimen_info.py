# -*- coding: utf-8 -*-
# __author__ = 'haidong.gu'

"""宏基因输入的specimen_info表，要求文件四列，有表头"""

from biocluster.iofile import File
from biocluster.core.exceptions import OptionError
from biocluster.core.exceptions import FileError
import re


class SpecimenInfoFile(File):
    """
    txt类
    """

    def __init__(self):
        super(SpecimenInfoFile, self).__init__()

    def check(self):
        """
        功能：实现对样本信息检查，如果样本中含有空行，去掉并返回
        :return:
        """
        f = open(self.prop['path'], 'r')
        lines = f.readlines()
        r = open(self.prop['path'], 'w')
        n = 0
        for line in lines:
            n += 1
            if n == 1:
                r.write(line.strip()+"\n")
            else:
                if line.strip() == "":
                    continue
                line_split = line.split("\t")
                if len(line_split) != 4:
                    raise FileError('文件应为四列，请检查', code="42800801")
                #if not re.search(r'[A-Z]+', line_split[0], re.I):
                    #raise FileError('样本名应包含字母', code="42800802")
                if re.search(r'[ ,-,.]', line_split[0], re.I):
                    raise FileError("样本名不可含有空格、中划线、句号，支持字母、数字、下划线", code="42800803")
                if not re.search(r'[A-Z]+', line_split[1], re.I):
                    raise FileError("样品来源为字符串，需含有字母", code="42800804")
                if re.search(r'[A-Z]+', line_split[2], re.I):
                    raise FileError("插入序列信息需要为数字，不可含字母", code="42800805")
                if re.search(r'[A-Z]+', line_split[3], re.I):
                    raise FileError("reads长度信息需要为数字，不可含字母", code="42800806")
                r.write(line.strip()+"\n")
        f.close()
        r.close()
        return True

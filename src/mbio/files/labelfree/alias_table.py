# -*- coding: utf-8 -*-
# __author__ = 'shicaiping'
import re
from biocluster.iofile import File
from biocluster.core.exceptions import FileError
from collections import defaultdict


class AliasTableFile(File):
    """
    定义alias_table格式文件
    """
    def __init__(self):
        super(AliasTableFile, self).__init__()

    def get_info(self):
        """
        获取文件属性
        """
        super(AliasTableFile, self).get_info()
        info = self.get_file_info()
        self.set_property("sample_number", len(info[0]))
        self.set_property("sample", info[0])

    def get_file_info(self):
        """
        获取alias_table文件的信息
        """
        row = 0
        self.format_check()
        with open(self.prop['path'], 'r') as f:
            sample = list()
            line = f.readline().rstrip()  # 将rstrip("\r\n") 全部替换为rstrip()
            line = re.split("\t", line)
            for line in f:
                line = line.rstrip()
                line = re.split("\t", line)
                row += 1
                if line[0] not in sample:
                    sample.append(line[0])
            return (sample)

    def format_check(self):
        with open(self.prop['path'], 'r') as f:
            line = f.readline().rstrip()
            if not re.search("^#", line[0]):
                raise FileError("该alias文件不含表头，alias表第一列应该以#号开头")
            line = line.split("\t")
            length = len(line)
            if length < 2:
                raise FileError('alias_table文件至少应该有两列')
            for i in line[1:]:
                if re.search("\s", i):
                    raise FileError('别名表里不可以包含空格')
        with open(self.prop['path'], 'r') as f:
            for line in f:
                if line.startswith("#"):
                    continue
                line = line.rstrip()
                line = re.split("\t", line)
                for l in line:
                    if re.search("\s", l):
                        raise FileError('别名表里不可以包含空格')
                len_ = len(line)
                if len_ != length:
                    raise FileError("文件的列数不相等")

    def check(self):
        if super(AliasTableFile, self).check():
            self.get_info()
            if self.prop['sample_number'] == 0:
                raise FileError('应该至少包含一个样本')

if __name__ == "__main__":
    g = AliasTableFile()
    g.set_path("example.alias")
    g.get_info()

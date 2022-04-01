# -*- coding: utf-8 -*-
# __author__ = 'zouxuan'
import re
from biocluster.iofile import File
from biocluster.core.exceptions import FileError


class InsertsizeTableFile(File):
    """
    定义插入片段文件格式
    """

    def __init__(self):
        super(InsertsizeTableFile, self).__init__()

    def get_info(self):
        """
        获取文件属性
        :return:
        """
        super(InsertsizeTableFile, self).get_info()
        info = self.get_file_info()
        self.set_property("sample_number", len(info[0]))
        self.set_property("sample", info[0])
        self.set_property("insertsize", info[1])

    def get_file_info(self):
        """
        获取insertsize_table文件的信息
        """
        row = 0
        self.format_check()
        with open(self.prop['path'], 'r') as f:
            sample = list()
            lenth = {}
            for line in f:
                line = line.rstrip()
                line = re.split("\t", line)
                row += 1
                if line[0] in sample:
                    raise FileError('样品名有重复', code="43900401")
                sample.append(line[0])
                lenth[line[0]] = line[1]
            return sample, lenth

    def format_check(self):
        with open(self.prop['path'], 'r') as f:
            for line in f:
                line = line.rstrip()
                line = re.split("\t", line)
                if isinstance(line[1], int):
                    raise FileError('插入片段长度必须为整数', code="43900402")
                if line[1] < 0:
                    raise FileError('插入片段长度必须大于0', code="43900403")
                for l in line:
                    if re.search("\s", l):
                        raise FileError('文件里不可以包含空格', code="43900404")
                len_ = len(line)
                if len_ != 2:
                    raise FileError('insertsize_table 文件应该为两列', code="43900405")

    def check(self):
        if super(InsertsizeTableFile, self).check():
            self.get_info()
            if self.prop['sample_number'] == 0:
                raise FileError('应该至少包含一个样本', code="43900406")


if __name__ == "__main__":
    i = InsertsizeTableFile()
    i.set_path("/mnt/ilustre/users/sanger-dev/sg-users/zouxuan/insertSize")

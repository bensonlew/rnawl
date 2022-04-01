# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'@20200617

from biocluster.iofile import File
from biocluster.core.exceptions import FileError

class RawSequenceTxtFile(File):
    """
    定义gff格式文件
    """

    def __init__(self):
        super(RawSequenceTxtFile, self).__init__()
        self.amplified_region = []
        self.sample_list = []
        self.insert_size = {}
        self.sequencing_length = {}
        self.raw_reads = {}
        self.total_base = {}

    def check(self):
        """
        检测文件是否满足要求，发生错误时应该触发FileError异常
        :return:
        """
        if super(RawSequenceTxtFile, self).check():
            info = self.get_info()
            self.set_property("column_number", info[0])
            self.set_property("sample_list", info[1])

    def get_info(self):
        super(RawSequenceTxtFile,self).get_info()
        with open(self.prop["path"],"r") as f:
            line = f.readline()
            column_number = 0
            if line.find("#") != 0:
                raise FileError("首行应该以#号开始")
            line = line.strip().split("\t")
            first_name = line[0]
            for line in f:
                line = line.strip()
                lst = line.split("\t")
                column_number = len(lst)
                if len(lst) != 6 and len(lst) != 5:
                    raise FileError("文件格式错误，原始序列文件应有5列或者6列")
                else:
                    if first_name in ["#amplified_region"]:
                        if not lst[0] in self.amplified_region:
                            self.amplified_region.append(lst[0])
                        else:
                            raise FileError("扩增子重复，请检查文件首列")
                    else:
                        if len(lst) == 6:
                            if not lst[0] in self.sample_list:
                                self.sample_list.append(lst[0])
                            if not lst[1] in self.amplified_region:
                                self.amplified_region.append(lst[1])
                        elif len(lst) == 5:
                            if not lst[1] in self.amplified_region:
                                self.amplified_region.append(lst[1])
                            self.sample_list.append("-")
                    self.insert_size[lst[0]] = lst[-4]
                    self.sequencing_length[lst[0]] = lst[-3]
                    self.raw_reads[lst[0]] = lst[-2]
                    self.total_base[lst[0]] = lst[-1]
        with open(self.prop["path"],"r") as m:
            line = m.readline()
            print(line)
            lines = m.readlines()
            if len(lines) <= 0:
                raise FileError("this file has only header,please check")
        return (column_number, self.sample_list)
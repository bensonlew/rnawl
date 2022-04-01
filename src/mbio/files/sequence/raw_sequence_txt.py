# -*- coding: utf-8 -*-
# __author__ = 'sj'

from biocluster.iofile import File
from biocluster.core.exceptions import FileError
import subprocess
import os

class RawSequenceTxtFile(File):
    """
    定义gff格式文件
    """

    def __init__(self):
        super(RawSequenceTxtFile, self).__init__()
        self.amplified_region = []
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
            self.get_info()
            
            
    def get_info(self):
        super(RawSequenceTxtFile,self).get_info()
        with open(self.prop["path"],"r") as f:
            line = f.readline()
            if line.find("#") != 0:
                raise FileError("首行应该以#号开始", code="44001601")
            for line in f:
                line = line.strip()
                lst = line.split("\t")
                # print lst
                if len(lst) != 5:
                    raise FileError("文件格式错误，原始序列文件应有5列", code="44001602")
                else:
                    if not lst[0] in self.amplified_region:
                        self.amplified_region.append(lst[0])
                    else:
                        raise FileError("扩增子重复，请检查文件首列", code="44001603")
                    self.insert_size[lst[0]] = lst[1]
                    self.sequencing_length[lst[0]] = lst[2]
                    self.raw_reads[lst[0]] = lst[3]
                    self.total_base[lst[0]] = lst[4]
        with open(self.prop["path"],"r") as m:
            line = m.readline()
            print(line)
            lines = m.readlines()
            if len(lines) <= 0:
                raise FileError("this file has only header,please check", code="44001604")
        return True

            
if __name__ == '__main__':
    a = RawSequenceTxtFile()
    a.set_path("/mnt/ilustre/users/sanger-dev/sg-users/shijin/raw_sequence.txt")
    a.check()
    print a.insert_size
    

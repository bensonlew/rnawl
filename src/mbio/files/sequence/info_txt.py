# -*- coding: utf-8 -*-
# __author__ = 'sj'

from biocluster.iofile import File
from biocluster.core.exceptions import FileError
import subprocess
import os

class InfoTxtFile(File):
    """
    定义gff格式文件
    """

    def __init__(self):
        super(InfoTxtFile, self).__init__()
        # self.workdir_path = []
        self.sample_path = []
        self.length_path = []

    def check(self):
        """
        检测文件是否满足要求，发生错误时应该触发FileError异常
        :return:
        """
        if super(InfoTxtFile, self).check():
            self.get_info()
            return True
            
    def get_info(self):
        super(InfoTxtFile,self).get_info()
        with open(self.prop["path"],"r") as f:
            line = f.readline()
            if line.find("#") != 0:
                raise FileError("首行应该以#号开始", code="44001201")
            for line in f:
                line = line.strip()
                lst = line.split("\t")
                # print lst
                if len(lst) != 8:
                    raise FileError("文件格式错误，原始序列文件应有5列", code="44001202")
                else:
                    # path = line.split("\t")[2]
                    path = lst[2]
                    # if path not in self.workdir_path:
                    #     self.workdir_path.append(path)
                    sample = path + "/output/fa/" + lst[1] + ".fasta"
                    length = path + "/output/length/" + lst[1] + ".length_file"
                    self.sample_path.append(sample)
                    self.length_path.append(length)
        return True

            

    

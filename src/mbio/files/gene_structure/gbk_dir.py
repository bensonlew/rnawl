# -*- coding: utf-8 -*-
# __author__ = gaohao
# last_modify：2018.03.16


import re
import os
import subprocess
from Bio import SeqIO
from biocluster.core.exceptions import FileError
from biocluster.iofile import Directory
from mbio.files.gene_structure.gbk import GbkFile

'''
微生物基因组GBK文件夹检查
'''


class GbkDirFile(Directory):
    """
    定义gbk_dir文件夹
    """
    def __init__(self):
        super(GbkDirFile, self).__init__()
        self.gbks =[]

    def get_info(self):
        """
        获取文件属性
        :return:
        """
        super(GbkDirFile, self).get_info()
        self.get_gbk_info()
        self.set_property("gbk_list", self.gbks)

    def check(self):
        """
        检测文件是否满足要求,发生错误时应该触发FileError异常
        :return:
        """
        filelist = os.listdir(self.prop['path'])
        if not len(filelist):
            raise FileError('原始序列文件夹为空，请检查确认', code="42200401")
        for file in filelist:
            if not re.search(r'.gbk',file):
                raise FileError('原始序列文件夹的文件不是以.gbk结尾，请检查确认', code="42200402")
            else:
                my_gbk = GbkFile()
                fq_path = os.path.join(self.prop['path'], file)
                my_gbk.set_path(fq_path)
                my_gbk.get_info()

    def get_gbk_info(self):
        """
        获取dir文件夹的样品信息
        :return: (sample_list)
        """
        filelist2 = os.listdir(self.prop['path'])
        print(self.prop['path'])
        for file in filelist2:
            print(file)
            my_gbk = GbkFile()
            fq_path = os.path.join(self.prop['path'], file)
            my_gbk.set_path(fq_path)
            my_gbk.get_info()
            if my_gbk.check():
                self.gbks.append(file)

if __name__ == "__main__":
    gbk_dir = GbkDirFile()

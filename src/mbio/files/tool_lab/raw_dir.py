# -*- coding: utf-8 -*-
# __author__ = gaohao
# last_modify：2018.03.16


import re, Bio, urllib2, regex, os
import subprocess
from biocluster.iofile import Directory
from collections import defaultdict
from biocluster.config import Config
from biocluster.core.exceptions import FileError
from mbio.files.sequence.fastq import FastqFile

'''
aimhii参数in_fastq文件检查
'''


class RawDirFile(Directory):
    """
    定义raw_dir文件夹
    """
    def __init__(self):
        super(RawDirFile, self).__init__()

    def get_info(self):
        """
        获取文件属性
        :return:
        """
        super(RawDirFile, self).get_info()
        seqinfo = self.get_raw_info()
        self.set_property("sample_list", seqinfo[0])

    def check(self):
        """
        检测文件是否满足要求,发生错误时应该触发FileError异常
        :return:
        """
        filelist = os.listdir(self.prop['path'])
        raw_list_path = os.path.join(self.prop['path'], "list.txt")
        sample_list = {}
        if not len(filelist):
            raise FileError('原始序列文件夹为空，请检查确认', code="41400801")
        if not os.path.exists(raw_list_path):
            raise FileError('原始序列文件夹的list.txt为不存在，请检查确认', code="41400802")
        else:
            with open(raw_list_path, "rb") as l:
                lines = l.readlines()
                for line in lines[1:]:  ##有表头
                    line2 = line.strip().split("\t")
                    if len(line2) == 2:
                        sample_list[line2[0]] = line2[1]
                        raw_path = line2[1].split(';')
                        for raw in raw_path:
                            files_path = raw.split(',')
                            for file_path in files_path:
                                if file_path.endswith(".fastq.gz"):
                                    if file_path not in filelist:
                                        raise FileError('原始序列文件夹的序列文件名称与list.txt不同，请检查确认')
                                else:
                                    raise FileError('测序数据必须为gunzip压缩文件')
                    else:
                        raise FileError('list.txt文件格式有误')
                if len(sample_list.values()) != len(set(sample_list.values())):
                    raise FileError('list.txt文件中有重复的序列文件名称，请检查确认')

    def get_raw_info(self):
        """
        获取dir文件夹的样品信息
        :return: (sample_list)
        """
        sample_list ={}
        raw_list_path = os.path.join(self.prop['path'], "list.txt")
        with open(raw_list_path, "rb") as l:
            lines = l.readlines()
            for line in lines[1:]:
                line2 = line.strip().split("\t")
                #if len(line2) == 6:
                sample_list[line2[0]] = line2[0]  ##
        return sample_list



if __name__ == "__main__":
    raw_dir = RawDirFile()
    #raw_dir.set_path("/mnt/ilustre/users/sanger-dev/sg-users/zhujuan/Genomic/database/")
    #test = gbk.get_info()
    #gbk.check()


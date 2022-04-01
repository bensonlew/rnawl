# -*- coding: utf-8 -*-
# __author__ = gaohao
# last_modify：2018.03.21


import re, Bio, urllib2, regex, os
import subprocess
from biocluster.iofile import Directory
from collections import defaultdict
from biocluster.config import Config
from biocluster.core.exceptions import FileError


'''
微生物基因组参数rawdata文件检查
'''


class SampleDirFile(Directory):
    """
    定义raw_dir文件夹
    """
    def __init__(self):
        super(SampleDirFile, self).__init__()

    def get_info(self):
        """
        获取文件属性
        :return:
        """
        super(SampleDirFile, self).get_info()
        seqinfo = self.get_raw_info()
        self.set_property("sample_list", seqinfo[0])

    def check(self):
        """
        检测文件是否满足要求,发生错误时应该触发FileError异常
        :return:
        """
        filelist = os.listdir(self.prop['path'])
        sample_list = {}
        if not len(filelist):
            raise FileError('原始序列文件夹为空，请检查确认', code="41400901")
        if "raw_list.txt" in filelist:
            raw_list_path = os.path.join(self.prop['path'], "raw_list.txt")
            with open(raw_list_path, "rb") as l:
                lines = l.readlines()
                for line in lines[1:]:
                    line2 = line.strip().split("\t")
                    if len(line2) == 6:
                        sample_list[line2[0]] = line2[0]
                        if re.search(';', line2[1]):
                            raw_path = line2[1].split(';')
                            for raw in raw_path:
                                files_path = raw.split(',')
                                for file_path in files_path:
                                    if file_path not in filelist:
                                        raise FileError('原始序列文件夹的序列名称与list.txt不同，请检查确认', code="41400902")
                        else:
                            files_path = line2[1].split(',')
                            for file_path in files_path:
                                if file_path not in filelist:
                                    raise FileError('原始序列文件夹的序列名称与list.txt不同，请检查确认', code="41400903")
                    else:
                        raise FileError('list.txt文件格式有误', code="41400904")

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
                if len(line2) == 6:
                    sample_list[line2[0]] = line2[0]
        return sample_list



if __name__ == "__main__":
    raw_dir = SampleDirFile()
    #raw_dir.set_path("/mnt/ilustre/users/sanger-dev/sg-users/zhujuan/Genomic/database/")
    #test = gbk.get_info()
    #gbk.check()

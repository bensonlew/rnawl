# -*- coding: utf-8 -*-
# __author__ = gaohao
# last_modify：2018.03.16


import re, Bio, urllib2, regex, os
import subprocess
from biocluster.iofile import Directory
from collections import defaultdict
from biocluster.config import Config
from biocluster.core.exceptions import FileError

'''
微生物基因组参数assemble文件夹检查
'''


class AsseDirFile(Directory):
    """
    定义raw_dir文件夹
    """
    def __init__(self):
        super(AsseDirFile, self).__init__()

    def get_info(self):
        """
        获取文件属性
        :return:
        """
        super(AsseDirFile, self).get_info()
        seqinfo = self.get_raw_info()
        self.set_property("sample_list", seqinfo[0])

    def check(self):
        """
        检测文件是否满足要求,发生错误时应该触发FileError异常
        :return:
        """
        filelist = os.listdir(self.prop['path'])
        asse_list_path = os.path.join(self.prop['path'], "list.txt")
        if not len(filelist):
            raise FileError('组装序列文件夹为空，请检查确认', code="42100101")
        if not os.path.exists(asse_list_path):
            raise FileError('组装序列文件夹的list.txt为不存在，请检查确认', code="42100102")


    def check_file(self):
        """
        当扫描图时，有rawdata时且上传组装assemble文件，检测文件是否满足要求,发生错误时应该触发FileError异常
        :return:
        """
        filelist = os.listdir(self.prop['path'])
        asse_list_path = os.path.join(self.prop['path'], "list.txt")
        with open(asse_list_path, "rb") as l:
            lines = l.readlines()
            for line in lines[1:]:
                line2 = line.strip('\r\n').split("\t")
                if len(line2) == 2:
                    if line2[1] not in filelist:
                        raise FileError('组装序列文件夹的序列文件名称与list.txt中的不同，请检查确认', code="42100103")
                else:
                    raise FileError('list.txt文件格式有误', code="42100104")

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
                sample_list[line2[0]] = line2[0]
        return sample_list



if __name__ == "__main__":
    raw_dir = AsseDirFile()
    #raw_dir.set_path("/mnt/ilustre/users/sanger-dev/sg-users/zhujuan/Genomic/database/")
    #test = file.get_info()
    #file.check()

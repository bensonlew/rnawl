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
微生物基因组参数assemble文件夹检查
'''


class FnaDirFile(Directory):
    """
    定义raw_dir文件夹
    """
    def __init__(self):
        super(FnaDirFile, self).__init__()

    def get_info(self):
        """
        获取文件属性
        :return:
        """
        super(FnaDirFile, self).get_info()
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
            raise FileError('组装序列文件夹为空，请检查确认', code="41400101")
        if not os.path.exists(asse_list_path):
            raise FileError('组装序列文件夹的list.txt为不存在，请检查确认', code="41400102")
        os.system('dos2unix %s'%asse_list_path)


    def check_file1(self):
        """
        当扫描图时，有rawdata时且上传组装assemble文件，检测文件是否满足要求,发生错误时应该触发FileError异常
        :return:
        """
        filelist = os.listdir(self.prop['path'])
        asse_list_path = os.path.join(self.prop['path'], "list.txt")
        with open(asse_list_path, "rb") as l:
            lines = l.readlines()
            for line in lines[1:]:
                line2 = line.strip().split("\t")
                if len(line2) == 2:
                    if line2[1] not in filelist:
                        raise FileError('组装序列文件夹的序列文件名称与list.txt中的不同，请检查确认', code="41400103")
                else:
                    raise FileError('list.txt文件格式有误', code="41400104")

    def check_file2(self):
        """
        当扫描图时，只有上传组装assemble文件，检测文件是否满足要求,发生错误时应该触发FileError异常
        :return:
        """
        filelist = os.listdir(self.prop['path'])
        asse_list_path = os.path.join(self.prop['path'], "list.txt")
        with open(asse_list_path, "rb") as l:
            lines = l.readlines()
            for line in lines[1:]:
                line2 = line.strip().split("\t")
                if len(line2) == 4:
                    if line2[3] not in filelist:
                        raise FileError('组装序列文件夹的序列文件名称与list.txt不同，请检查确认', code="41400105")
                else:
                    raise FileError('list.txt文件格式有误', code="41400106")

    def check_file3(self):
        """
        当完成图时，有rawdata时且上传组装assemble文件，检测文件是否满足要求,发生错误时应该触发FileError异常
        :return:
        """
        filelist = os.listdir(self.prop['path'])
        asse_list_path = os.path.join(self.prop['path'], "list.txt")
        with open(asse_list_path, "rb") as l:
            lines = l.readlines()
            for line in lines[1:]:
                line2 = line.strip().split("\t")
                if len(line2) == 3:
                    if line2[1] not in filelist:
                        raise FileError('组装序列文件夹的序列文件名称与list.txt不同，请检查确认', code="41400107")
                    if line2[2] not in ['chromosome','Chromosome','plasmid','Plasmid']:
                        raise FileError('基因组类型不正确！', code="41400108")
                else:
                    raise FileError('list.txt文件格式有误', code="41400109")

    def check_file4(self):
        """
        当完成图时，有rawdata时且上传组装assemble文件，检测文件是否满足要求,发生错误时应该触发FileError异常
        :return:
        """
        filelist = os.listdir(self.prop['path'])
        asse_list_path = os.path.join(self.prop['path'], "list.txt")
        with open(asse_list_path, "rb") as l:
            lines = l.readlines()
            for line in lines[1:]:
                line2 = line.strip().split("\t")
                if len(line2) == 5:
                    if line2[3] not in filelist:
                        raise FileError('组装序列文件夹的序列文件名称与list.txt不同，请检查确认', code="41400110")
                    if line2[4] not in ['chromosome','Chromosome','plasmid','Plasmid']:
                        raise FileError('基因组类型不正确！', code="41400111")
                else:
                    raise FileError('list.txt文件格式有误', code="41400112")

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
    raw_dir = FnaDirFile()
    #raw_dir.set_path("/mnt/ilustre/users/sanger-dev/sg-users/zhujuan/Genomic/database/")
    #test = file.get_info()
    #file.check()

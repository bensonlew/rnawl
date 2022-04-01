# -*- coding: utf-8 -*-
# __author__ = gaohao
# last_modify：2018.06.03


import re, Bio, urllib2, regex, os
import subprocess
from biocluster.iofile import Directory
from collections import defaultdict
from biocluster.config import Config
from biocluster.core.exceptions import FileError
from mbio.files.sequence.fastq import FastqFile

'''
微生物基因组pacbio_dir文件夹
'''


class PacbioDirFile(Directory):
    """
    定义raw_dir文件夹
    """
    def __init__(self):
        super(PacbioDirFile, self).__init__()

    def get_info(self):
        """
        获取文件属性
        :return:
        """
        super(PacbioDirFile, self).get_info()
        seqinfo = self.get_raw_info()
        self.set_property("sample_list", seqinfo[0])

    def check(self):
        """
        检测文件是否满足要求,发生错误时应该触发FileError异常
        :return:
        """
        filelist = os.listdir(self.prop['path'])
        raw_list_path = os.path.join(self.prop['path'], "pacbio.rawdata.list")
        sample_list = {}
        if not len(filelist):
            raise FileError('原始序列文件夹为空，请检查确认', code="41400701")
        if not os.path.exists(raw_list_path):
            raise FileError('原始序列文件夹的pacbio.rawdata.list为不存在，请检查确认', code="41400702")
        else:
            with open(raw_list_path, "rb") as l:
                lines = l.readlines()
                for line in lines[1:]:
                    line2 = line.strip().split("\t")
                    if len(line2) == 6:
                        sample_list[line2[0]] = line2[0]
                        if line2[1] not in filelist:
                            raise FileError('文件夹的序列文件名称与pacbio.rawdata.list不同，请检查确认', code="41400703")
                    else:
                        raise FileError('list.txt文件格式有误', code="41400704")
                    if line2[5] not in ['pacbio','Pacbio']:
                        raise FileError('list.txt文件的lib是有错的！', code="41400705")
                    if re.search(r'M$',line2[4]):
                        raise FileError('请提供基因组大小且单位为M！', code="41400706")

    def get_raw_info(self):
        """
        获取dir文件夹的样品信息
        :return: (sample_list)
        """
        sample_list ={}
        raw_list_path = os.path.join(self.prop['path'], "pacbio.rawdata.list")
        with open(raw_list_path, "rb") as l:
            lines = l.readlines()
            for line in lines[1:]:
                line2 = line.strip().split("\t")
                if len(line2) == 6:
                    sample_list[line2[0]] = line2[0]
        return sample_list



if __name__ == "__main__":
    raw_dir = PacbioDirFile()
    #raw_dir.set_path("/mnt/ilustre/users/sanger-dev/sg-users/zhujuan/Genomic/database/")
    #test = gbk.get_info()
    #gbk.check()


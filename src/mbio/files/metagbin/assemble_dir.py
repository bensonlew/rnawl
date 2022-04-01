# -*- coding: utf-8 -*-
# __author__ = 'xuting'
import re
import os
from biocluster.core.exceptions import FileError
from biocluster.iofile import Directory
import operator

class AssembleDirFile(Directory):
    """
    定义fasta文件夹
    需要biopython
    """
    def __init__(self):
        """
        """
        super(AssembleDirFile, self).__init__()

    def check_assmble(self):
        filelist = os.listdir(self.prop['path'])
        files=[]
        files2=[]
        for file in filelist:
            if not file.endswith((".gz", ".txt", ".tar.gz", "stat.xls")):
                raise FileError('组装序列文件夹的文件有误')
            else:
                if file.endswith((".gz", ".tar.gz")):
                    files2.append(file)
        with open (os.path.join(self.prop['path'], "list.txt"),'r') as f:
            lines =f.readlines()
            for lin in lines[1:]:
                line =lin.rstrip('\r\n').split('\t')
                files.append(line[1])
        files = sorted(files)
        files2 = sorted(files2)
        if not operator.eq(files,files2):
            raise FileError('组装文件夹文件与list文件的文件名称不一样！')

    def check(self):
        """
        检测文件夹是否满足要求，不满足时触发FileError异常
        :return:
        """
        if super(AssembleDirFile, self).check():
            self.check_assmble()
            return True
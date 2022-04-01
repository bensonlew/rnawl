# -*- coding: utf-8 -*-
# __author__ = 'yuguo'
# creat at 20171111

from biocluster.iofile import Directory
from biocluster.core.exceptions import FileError
import os
import re


class SampleDirFile(Directory):
    '''
    多个样本fastq.gz文件目录
    '''
    def __init__(self):
        super(SampleDirFile, self).__init__()

    def get_info(self):
        '''
        获取文件属性
        '''
        super(SampleDirFile, self).get_info()

    def check(self):
        '''
        检查文件格式
        '''
        if super(SampleDirFile, self).check():
            self.get_info()
            flist = os.listdir(self.prop['path'])
            if len(flist) == 0:
                raise FileError("文件夹为空")
            # for f in flist:
            #     m = re.match('.*(fq|fastq)(\.gz|)', f)
            #     if not m:
            #         raise FileError("样本文件格式不正确")
        else:
            raise FileError("文件格式错误")
        return True

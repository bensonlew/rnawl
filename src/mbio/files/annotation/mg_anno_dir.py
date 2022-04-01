# -*- coding: utf-8 -*-
# __author__ = 'shaohua.yuan'
import re
import os
import subprocess
from biocluster.config import Config
from biocluster.iofile import Directory
from biocluster.core.exceptions import FileError
from mbio.files.sequence.fastq import FastqFile
from mbio.files.sequence.file_sample import FileSampleFile


class MgAnnoDirFile(Directory):
    """
    定义宏基因nr注释结果文件夹的格式
    """
    def __init__(self):
        super(MgAnnoDirFile, self).__init__()

    def check_info(self):
        """
        获取判断文件夹内容
        """
        if 'path' in self.prop.keys() and os.path.isdir(self.prop['path']):
            # file_name = os.listdir(self.prop["path"])
            # if len(file_name) != 8:
            #     raise FileError("文件夹中内容不对,请检查!")
            # else:
            pass
        else:
            raise FileError("文件夹路径不正确，请设置正确的文件夹路径!", code="41200901")

    def check(self):
        """
        检测文件夹是否满足要求，不满足时触发FileError异常
        :return:
        """
        if super(MgAnnoDirFile, self).check():
            self.check_info()
            return True

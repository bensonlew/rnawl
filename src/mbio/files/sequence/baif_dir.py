# -*- coding: utf-8 -*-
# __author__ = 'zhouxuan'
import re
import os
import subprocess
from biocluster.config import Config
from biocluster.iofile import Directory
from biocluster.core.exceptions import FileError
from mbio.files.sequence.fastq import FastqFile
from mbio.files.sequence.file_sample import FileSampleFile


class BaifDirFile(Directory):
    """
    定义base_info文件夹的格式
    """
    def __init__(self):
        super(BaifDirFile, self).__init__()

    def get_info(self):
        """
        获取判断文件夹内容
        """
        if 'path' in self.prop.keys() and os.path.isdir(self.prop['path']):
            file_name = os.listdir(self.prop["path"])
            if len(file_name) == 0:
                raise FileError("文件夹为空，请设置正确的文件夹路径！", code="44000101")
            for name in file_name:
                file_path = os.path.join(self.prop['path'], name)
                with open(file_path, 'r') as r:
                    for line in r:
                        if len(line.strip("\n").split("\t")) == 18:
                            pass
                        else:
                            raise FileError("文件夹中%s文件内容不正确，请确认！", variables=(name), code="4400102")
            # self.set_property("stat_number", self.get_fastq_number())
        else:
            raise FileError("文件夹路径不正确，请设置正确的文件夹路径!", code="44000103")

    def check(self):
        """
        检测文件夹是否满足要求，不满足时触发FileError异常
        :return:
        """
        if super(BaifDirFile, self).check():
            self.get_info()
            return True

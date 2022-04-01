# -*- coding: utf-8 -*-

from biocluster.iofile import File
from biocluster.core.exceptions import FileError
import subprocess
import os

class BamlistFile(File):
    """
    定义Bamlist格式文件
    """

    def __init__(self):
        super(BamlistFile, self).__init__()

    def get_info(self):
        super(BamlistFile, self).get_info()


    def check(self):
        if super(BamlistFile,self).check():
            self.get_info()
            return True


if __name__ == '__main__':
    a = BamlistFile()
    a.set_path("/mnt/ilustre/users/sanger-dev/sg-users/liubinxu/test_denovo/test_data2/bam.list")
    a.check()

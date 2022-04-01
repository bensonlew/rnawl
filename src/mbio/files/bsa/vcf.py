# -*- coding: utf-8 -*-

from biocluster.iofile import File
from biocluster.core.exceptions import FileError
import subprocess
import os
# 写一个空文件让他任意检查通过


class VcfFile(File):
    """
    定义Vcf格式文件
    """

    def __init__(self):
        super(VcfFile, self).__init__()

    def get_info(self):
        super(VcfFile, self).get_info()

    def check(self):
        if super(VcfFile,self).check():
            self.get_info()
            return True


if __name__ == '__main__':
    a = VcfFile()
    a.set_path("/mnt/ilustre/users/sanger-dev/sg-users/litangjian/SNP/call.vcf")
    a.check()

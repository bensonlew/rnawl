# -*- coding: utf-8 -*-
# __author__ = 'zhaobinbin'

import re
import os
from biocluster.iofile import Directory
from biocluster.core.exceptions import FileError


class WtDirFile(Directory):
    """
    检测选择性消除的结果文件夹
    """
    def __init__(self):
        super(WtDirFile, self).__init__()

    def is_file(self, file_path):
        """
        检查是否是文件是否存在
        :param file_path:
        :return:
        """
        if not os.path.isfile(file_path) or not os.path.exists(file_path):
            raise FileError("{}中不存在{}文件,请检查！".format(self.prop['path'], file_path))

    def is_dir(self, dir_path):
        """
        检查文件夹是否存在
        :param dir_path:
        :return:
        """
        if not os.path.isdir(dir_path) or not os.path.exists(dir_path):
            raise FileError("不存在{}路径,请检查！".format(dir_path))

    def check(self):
        """
        检测文件夹是否满足要求，不满足时触发FileError异常
        :return:
        """
        if super(WtDirFile, self).check():
            # bed = os.path.join(self.prop['path'], "pop.bed")
            # bim = os.path.join(self.prop['path'], "pop.bim")
            # fam = os.path.join(self.prop['path'], "pop.fam")
            # self.is_file(bed)
            # self.is_file(bim)
            # self.is_file(fam)
            return True

if __name__ == "__main__":
    a = WtDirFile()
    a.set_path("/mnt/ilustre/users/sanger-dev/sg-users/xuanhongdong/evolution/pop")
    print a.check()
    print "检查通过"

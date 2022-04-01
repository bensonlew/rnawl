# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# modified 2018.08.22

import os
from biocluster.iofile import File
from biocluster.core.exceptions import FileError


class PiTajimadFstFile(File):
    """
    定义pi_tajimaD_fst.selec格式
    必须是七列，不能为空文件
    """
    def __init__(self):
        super(PiTajimadFstFile, self).__init__()

    def check(self):
        if super(PiTajimadFstFile, self).check():
            with open(self.prop["path"], "r") as f:
                head = f.readline().strip().split(" ")
                if len(head) != 7:
                    raise FileError("{}必须是七列".format(self.prop["path"]))
                try:
                    line = f.next()
                except:
                    raise FileError("{}文件内容为空，请检查".format(self.prop["path"]))
            return True


if __name__ == "__main__":
    a = PiTajimadFstFile()
    a.set_path("/mnt/ilustre/users/sanger-dev/sg-users/zengjing/dna_evolution/file/4-6.1.pi_tajimaD_fst.select")
    a.check()

# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# modified 2018.08.22

import os
from biocluster.iofile import File
from biocluster.core.exceptions import FileError


class WindowPiFile(File):
    """
    定义windowed.pi格式
    必须是五列，不能为空文件
    """
    def __init__(self):
        super(WindowPiFile, self).__init__()

    def check(self):
        if super(WindowPiFile, self).check():
            with open(self.prop["path"], "r") as f:
                head = f.readline().strip().split("\t")
                if len(head) != 5:
                    raise FileError("{}必须是5列".format(self.prop["path"]))
                try:
                    line = f.next()
                    self.set_property("null", False)
                except:
                    self.set_property("null", True)
                    raise FileError("{}文件内容为空，请检查".format(self.prop["path"]))
            return True


if __name__ == "__main__":
    a = WindowPiFile()
    a.set_path("/mnt/ilustre/users/sanger-dev/sg-users/zengjing/dna_evolution/sweep/1..windowed.pi")
    a.check()

# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# modified 2018.08.22

import os
from biocluster.iofile import File
from biocluster.core.exceptions import FileError


class WindowWeirFstFile(File):
    """
    定义windowed.weir.fst格式
    必须是六列，不能为空文件
    """
    def __init__(self):
        super(WindowWeirFstFile, self).__init__()

    def check(self):
        if super(WindowWeirFstFile, self).check():
            with open(self.prop["path"], "r") as f:
                head = f.readline().strip().split("\t")
                if len(head) != 6:
                    raise FileError("{}必须是六列".format(self.prop["path"]))
                try:
                    line = f.next()
                    self.set_property("null", False)
                except:
                    self.set_property("null", True)
                    print "{}文件内容为空".format(self.prop["path"])
            return True


if __name__ == "__main__":
    a = WindowWeirFstFile()
    a.set_path("/mnt/ilustre/users/sanger-dev/workspace/20180822/Single_vcftools_stat2/VcftoolsStat/output/1-2.windowed.weir.fst")
    a.check()

# -*- coding: utf-8 -*-
# __author__: zengjing
# modified: 20190225

import os
from biocluster.iofile import File
from biocluster.core.exceptions import FileError


class VcfFile(File):
    """
    vcf文件检查,vcf文件不能为空
    """
    def __init__(self):
        super(VcfFile, self).__init__()

    def check(self):
        super(VcfFile, self).check()
        if not os.path.exists(self.prop["path"]):
            raise FileError("文件：%s不存在，请检查" % self.prop["path"])
        is_null = True
        with open(self.prop["path"], "r") as f:
            while is_null:
                line = f.readline()
                if not line.startswith("#"):
                    try:
                        line = f.readline()
                        is_null = False
                    except:
                        raise FileError("文件：%s是空的，请检查" % self.prop["path"])


if __name__ == "__main__":
    a = VcfFile()
    a.set_path("/mnt/ilustre/users/isanger/sg-users/zengjing/ref_wgs_v2/sv/test.out.vcf")
    a.check()

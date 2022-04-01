# -*- coding: utf-8 -*-
# __author__ = 'wangzhaoyue'

"""针对非多样性barcode拆分格式文件类"""

from biocluster.iofile import File
from biocluster.core.exceptions import FileError


class BarcodeInfoFile(File):
    """
   针对非多样性barcode拆分格式文件类
    """
    def __init__(self):
        super(BarcodeInfoFile, self).__init__()

    def check(self):
        """
        检测文件是否满足要求，发生错误时应该触发FileError异常
        :return:
        """
        if super(BarcodeInfoFile, self).check():
            with open(self.prop["path"], "r") as f:
                lines = f.readlines()
                tmp_len = lines[0].strip().split("\t")
                for line in lines[1:]:
                    tmp = line.strip().split("\t")
                    if len(tmp) != len(tmp_len):
                        raise FileError("%s行文件信息不全，应有%s列，实际%s列，请核实", variables=(line, len(tmp_len), len(tmp)), code="44000301")
                    # if tmp[4] not in ["PE", "RAD", "GBS"]:
                    #     raise FileError("文库类型{}不在[PE,RAD,GBS]中，请核实".format(tmp[5]))
        else:
            raise FileError("文件格式错误", code="44000302")

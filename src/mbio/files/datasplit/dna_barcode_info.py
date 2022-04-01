# -*- coding: utf-8 -*-
# __author__ = "zengjing"
# last_modify: 20190124

import os
import re
from biocluster.iofile import File
from collections import defaultdict
from biocluster.core.exceptions import FileError


class DnaBarcodeInfoFile(File):
    """
    定义DNAbarcode信息表
    RAD：enzyme1序列\tsample
    GBS: enzyme1序列\tenzyme2序列\tsample
    """
    def __init__(self):
        super(DnaBarcodeInfoFile, self).__init__()

    def check(self):
        super(DnaBarcodeInfoFile, self).check()
        if not os.path.exists(self.prop["path"]):
            raise FileError("文件:{}不存在，请检查".format(self.prop["path"]))
        lib_type_list = []
        with open(self.prop["path"], "r") as f:
            lines = f.readlines()
            for line in lines:
                item = line.strip().split("\t")
                if len(item) == 2:
                    lib_type = "RAD"
                    for k in item[0]:
                        if k not in ["A", "T", "C", "G"]:
                            raise FileError("{}里有酶的序列不为ATCG碱基序列，请检查".format(self.prop["path"]))
                elif len(item) == 3:
                    lib_type = "GBS"
                    for k in item[0]:
                        if k not in ["A", "T", "C", "G"]:
                            raise FileError("{}里有酶的序列不为ATCG碱基序列，请检查".format(self.prop["path"]))
                    for k in item[1]:
                        if k not in ["A", "T", "C", "G"]:
                            raise FileError("{}里有酶的序列不为ATCG碱基序列，请检查".format(self.prop["path"]))
                else:
                    raise FileError("{}barcode信息错误，必须是两列或者三列".format(self.prop["path"]))
                lib_type_list.append(lib_type)
        if len(list(set(lib_type_list))) != 1:
            raise FileError("文件:{}的文库类型不唯一，请检查".format(self.prop["path"]))

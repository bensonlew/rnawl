# -*- coding: utf-8 -*-
# __author__ = "zengjing"
# last_modify: 20190128

import os
import re
from biocluster.iofile import File
from collections import defaultdict
from biocluster.core.exceptions import FileError


class MetaBarcodeInfoFile(File):
    """
    定义META barcode信息表
    #Sample\tBarcode-tag\tFbarcode\tRbarcode\n
    """
    def __init__(self):
        super(MetaBarcodeInfoFile, self).__init__()

    def check(self):
        super(MetaBarcodeInfoFile, self).check()
        if not os.path.exists(self.prop["path"]):
            raise FileError("文件:{}不存在，请检查".format(self.prop["path"]))
        sample_list, barcode_name_list, barcode_seq_list = [], [], []
        with open(self.prop["path"], "r") as f:
            lines = f.readlines()
            for line in lines:
                item = line.strip().split("\t")
                if item[0] in sample_list:
                    raise FileError("{}里的样本：{}重复，请检查".format(self.prop["path"], item[0]))
                sample_list.append(item[0])
                if item[1] in barcode_name_list:
                    raise FileError("{}里的barcode：{}重复，请检查".format(self.prop["path"], item[1]))
                barcode_name_list.append(item[1])
                if len(item) == 4:
                    seq = item[2] + "_" + item[3]
                elif len(item) == 3:
                    seq = item[2]
                else:
                    raise FileError("{}barcode信息错误，只能是三列或四列".format(self.prop["path"]))
                if seq in barcode_seq_list:
                    raise FileError("{}里的barcode：{}重复，请检查".format(self.prop["path"], seq))
                barcode_seq_list.append(seq)

# -*- coding: utf-8 -*-
# __author__ = "zengjing"
# last_modify: 20190128

import os
import re
from biocluster.iofile import File
from collections import defaultdict
from biocluster.core.exceptions import FileError


class MetaPrimerInfoFile(File):
    """
    定义META primer信息表
    #Sample\tF-barcode\tLinkPrimer\tR-barcode\tReversePrimer\n
    """
    def __init__(self):
        super(MetaPrimerInfoFile, self).__init__()

    def check(self):
        super(MetaPrimerInfoFile, self).check()
        if not os.path.exists(self.prop["path"]):
            raise FileError("文件:{}不存在，请检查".format(self.prop["path"]))
        sample_list, barcode_list = [], []
        with open(self.prop["path"], "r") as f:
            lines = f.readlines()
            for line in lines:
                item = line.strip().split("\t")
                if len(item) != 5:
                    raise FileError("{}barcode信息错误，必须是五列".format(self.prop["path"]))
                if item[0] in sample_list:
                    raise FileError("{}里的样本：{}重复，请检查".format(self.prop["path"], item[0]))
                sample_list.append(item[0])
                seq = item[1] + "_" + item[3]
                # if seq in barcode_list:
                #     raise FileError("{}里的barcode：{}重复，请检查".format(self.prop["path"], seq))
                barcode_list.append(seq)

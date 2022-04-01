# -*- coding: utf-8 -*-
# __author__ = "zengjing"
# last_modify: 20190128

import os
import re
from biocluster.iofile import File
from collections import defaultdict
from biocluster.core.exceptions import FileError


class MicroSampleInfoFile(File):
    """
    定义微生物基因组样本信息表
    #Sample\tLibrary\tLibrary_type\tInsertSize\tPath\n
    """
    def __init__(self):
        super(MicroSampleInfoFile, self).__init__()

    def check(self):
        super(MicroSampleInfoFile, self).check()
        if not os.path.exists(self.prop["path"]):
            raise FileError("文件:{}不存在，请检查".format(self.prop["path"]))
        sample_list = []
        with open(self.prop["path"], "r") as f:
            lines = f.readlines()
            for line in lines:
                item = line.strip().split("\t")
                name = item[1] + ":" + item[0]
                if len(item) < 5:
                    raise FileError("{}文件必须大于等于5列，请检查".format(self.prop["path"]))
                if name in sample_list:
                    raise FileError("{}中文库样本：{}信息重复，请检查".format(self.prop["path"], name))
                sample_list.append(name)

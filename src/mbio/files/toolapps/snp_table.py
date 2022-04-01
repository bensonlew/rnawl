# -*- coding: utf-8 -*-

import os
from biocluster.iofile import File
from biocluster.core.exceptions import FileError


class SnpTableFile(File):
    """
    Binbin Zhao@20180227
    """

    def __init__(self):
        super(SnpTableFile, self).__init__()

    def is_exists(self):
        if not os.path.isfile(self.path) or not os.path.exists(self.path):
            raise FileError("原始文件中不存在{}文件！".format(self.path))

    def check(self):
        specimen_ids = []
        group_info = {}
        if super(SnpTableFile, self).check():
            self.is_exists()
            # with open(self.path, "r") as f:
            #     lines = f.readlines()
            #     for line in lines:
            #         col = line.strip().split(":")
            #         if len(col) != 2:
            #             col = line.strip().split("\t")
            #             if len(col) != 2:
            #                 raise FileError("文件:{}必须为两列，请检查".format(self.path))
            #         if col[1] not in group_info.keys():
            #             group_info[col[1]] = []
            #             specimen_ids.append(col[0])
            #         group_info[col[1]].append(col[0])
            # if not group_info:
            #     raise FileError("没有得到具体的分组方案，请检查分组文件:{}是否为空".format(self.path))
            # self.set_property("specimen_ids", specimen_ids)
            # self.set_property("group_info", group_info)
            return True

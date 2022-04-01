#-*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'
# last_modified: 20190507
from biocluster.core.exceptions import FileError
from biocluster.iofile import Directory
import os, re

class JsonDirFile(Directory):
    def __init__(self):
        """
        json文件夹
        :return:
        """
        super(JsonDirFile, self).__init__()

    def check(self):
        """
        检测文件夹是否满足要求，不满足时触发FileError异常
        :return:
        """
        if super(JsonDirFile, self).check():
            self.get_info()
            if not os.path.exists(self.prop["path"]):
                raise FileError("{}文件不存在".format(self.prop["path"]))
            self.check_info()
            return True

    def get_info(self):
        """
        获取文件属性
        """
        super(JsonDirFile, self).get_info()

    def check_info(self):
        """
        检查json文件夹
        :return:
        """
        json_dir = self.prop['path']
        for file in os.listdir(json_dir):
            if not file.endswith("json"):
                raise FileError('文件不是json格式')
        return True

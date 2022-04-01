# -*- coding: utf-8 -*-
# __author__ = 'xuting'
import re
from biocluster.iofile import File
from biocluster.core.exceptions import FileError


class VennTableFile(File):
    """
    定义Venn_Table
    """
    def __init__(self):
        """
        :param group:分组的名称
        """
        super(VennTableFile, self).__init__()
        self.group = list()

    def get_info(self):
        """
        获取文件属性
        """
        super(VennTableFile, self).get_info()
        self.set_property("group", self.get_file_info())

    def get_file_info(self):
        """
        获取表格文件的信息
        """
        if self.is_set:
            with open(self.prop['prop'], 'r') as f:
                line = f.readline()
                line = re.split(r'\t', line)
                if re.search(r'(.+)\s+only', line[0]):
                    self.group.append(re.search(r'(.+)\s+only', line[0]).group(1))

    def check(self):
        """
        检测文件格式
        """
        if super(VennTableFile, self).check():
            self.get_info()
            with open(self.prop['path']) as f:
                line = f.readline().rstrip('\n')
                line = re.split('\t')
                if len(line) != 3:
                    raise FileError('文件格式错误')
            return True

# -*- coding: utf-8 -*-
# __author__ = 'xuting'
from biocluster.iofile import File
from biocluster.core.exceptions import FileError


class PanCoreTableFile(File):
    """
    定义pan_core Otu表
    """
    def __init__(self):
        super(PanCoreTableFile, self).__init__()

    def get_info(self):
        """
        获取文件属性
        """
        super(PanCoreTableFile, self).get_info()

    def check(self):
        """
        检测文件是否满足要求,发生错误时应该触发FileError异常
        """
        if super(PanCoreTableFile, self).check():
            self.get_info()
            self.check_format()
            return True

    def check_format(self):
        """
        检查文件格式
        """
        row = 0
        with open(self.prop['path'], 'r') as f:
            for line in f:
                row += 1
            if row < 2:
                raise FileError("文件格式错误")

if __name__ == "__main__":
    a = PanCoreTableFile()
    a.set_path("core.richness.xls")

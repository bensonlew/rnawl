# -*- coding: utf-8 -*-

"""txt格式文件类"""

from biocluster.iofile import File
from biocluster.core.exceptions import OptionError


class VennTableFile(File):
    def __init__(self):
        super(VennTableFile, self).__init__()
        
    def check(self):
        """
        检测文件是否满足要求，发生错误时应该触发FileError异常
        :return:
        """
        if super(VennTableFile, self).check():
            return True
        else:
            raise FileError("文件格式错误")

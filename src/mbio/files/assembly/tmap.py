# -*- coding: utf-8 -*-
# __author__ = 'wangzhaoyue'

"""txt格式文件类"""

from biocluster.iofile import File
from biocluster.core.exceptions import OptionError


class TmapFile(File):
    """
    tmap类
    """
    def __init__(self):
        super(TmapFile, self).__init__() 
     
    def check(self):
        """
        检测文件是否满足要求，发生错误时应该触发FileError异常
        :return:
        """
        if super(TmapFile, self).check():
            super(TmapFile,self).get_info()
            with open(self.prop["path"],"r") as f:
                for line in f:
                    if line.find("#") == -1:
                        line = line.strip()
                        lst = line.split("\t")
                        if len(lst) !=13:
                            raise FileError("文件格式错误，tmap应有13列")
                        else:
                            return True
        else:
            raise FileError("文件格式错误")

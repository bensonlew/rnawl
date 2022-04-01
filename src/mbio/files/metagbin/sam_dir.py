# -*- coding: utf-8 -*-
# __author__ = 'gaohao'

"""bam类"""

from biocluster.iofile import Directory
from biocluster.core.exceptions import FileError
import os

class SamDirFile(Directory):
    """
    sam文件夹
    """
    def __init__(self):
        super(SamDirFile, self).__init__()

    def get_info(self):
        """
        获取文件属性
        """
        super(SamDirFile, self).get_info()

    def check(self):
        """
        检测文件是否满足要求，发生错误时应该触发FileError异常
        """
        if super(SamDirFile, self).check():
            self.get_info()
            filelist = os.listdir(self.prop['path'])
            for file in filelist:
                if file.endswith('.sam'):
                    pass
                else:
                    raise FileError("文件夹中必须都是sam格式")
            return True
        else:
            raise FileError("文件格式错误")
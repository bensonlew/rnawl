# -*- coding: utf-8 -*-
# __author__ = 'shijin'

"""star_index类"""

from biocluster.iofile import Directory
from biocluster.core.exceptions import FileError

class StarIndexFile(Directory):
    """
    bam文件夹
    """
    def __init__(self):
        super(StarIndexFile, self).__init__()

    def get_info(self):
        """
        获取文件属性
        """
        super(StarIndexFile, self).get_info()

    def check(self):
        """
        检测文件是否满足要求，发生错误时应该触发FileError异常
        """
        if super(StarIndexFile, self).check():
            return True
        else:
            raise FileError("文件格式错误")



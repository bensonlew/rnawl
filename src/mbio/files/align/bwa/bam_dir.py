# -*- coding: utf-8 -*-
# __author__ = 'wangzhaoyue'

"""bam类"""

from biocluster.iofile import Directory
from biocluster.core.exceptions import FileError
import os

class BamDirFile(Directory):
    """
    bam文件夹
    """
    def __init__(self):
        super(BamDirFile, self).__init__()

    def get_info(self):
        """
        获取文件属性
        """
        super(BamDirFile, self).get_info()

    def check(self):
        """
        检测文件是否满足要求，发生错误时应该触发FileError异常
        """
        if super(BamDirFile, self).check():
            self.get_info()
            filelist = os.listdir(self.prop['path'])
            for file in filelist:
                if file.endswith('.bam'):
                    pass
                else:
                    raise FileError("文件夹中必须都是bam格式的文件", code = "41100801")
            return True
        else:
            raise FileError("文件格式错误", code = "41100802")
# if __name__ == "__main__":
#     a = BamDirFile()
#     a.set_path('input')
#     a.check()


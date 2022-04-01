# -*- coding: utf-8 -*-
# __author__ = 'HONGDONG'


from biocluster.iofile import Directory


class DirFile(Directory):
    """
    空的，只适用于对象存储下载文件用
    """
    def __init__(self):
        super(DirFile, self).__init__()

    def check(self):
        """
        检测文件夹是否满足要求，不满足时触发FileError异常
        :return:
        """
        if super(DirFile, self).check():
            return True

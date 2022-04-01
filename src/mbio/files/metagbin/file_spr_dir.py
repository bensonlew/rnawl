# -*- coding: utf-8 -*-
# __author__ = gaohao
# last_modify：2020.04.20


from biocluster.iofile import Directory
import os,re
from biocluster.core.exceptions import FileError

'''
spr的文件夹带list.txt准备文件夹的检查
'''


class FileSprDirFile(Directory):
    """
    """

    def __init__(self):
        super(FileSprDirFile, self).__init__()

    def get_info(self):
        """
        获取文件属性
        :return:
        """
        super(FileSprDirFile, self).get_info()
        file_list = os.listdir(self.prop['path'])
        for file in file_list:
            if not file.endswith((".spr", ".txt")):
                raise FileError('%s该目录下不是.spr压缩文件或list.txt，请检查！', variables=(file))

    def check(self):
        if super(FileSprDirFile, self).check():
            self.get_info()
            return True


if __name__ == "__main__":
    dir = FileSprDirFile()
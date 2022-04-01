# -*- coding: utf-8 -*-
# __author__ = gaohao
# last_modify：2019.01.15


from biocluster.iofile import Directory
import os,re
from biocluster.core.exceptions import FileError

'''
amphora准备文件夹的检查
'''


class AmphoraDirFile(Directory):
    """
    """

    def __init__(self):
        super(AmphoraDirFile, self).__init__()

    def get_info(self):
        """
        获取文件属性
        :return:
        """
        super(AmphoraDirFile, self).get_info()
        file_list = os.listdir(self.prop['path'])
        for file in  file_list:
            if not re.search('.anno.xls',file):
                raise FileError('%s文件不正确，请检查！', variables=(file))

    def check(self):
        if super(AmphoraDirFile, self).check():
            self.get_info()
            return True


if __name__ == "__main__":
    dir = AmphoraDirFile()
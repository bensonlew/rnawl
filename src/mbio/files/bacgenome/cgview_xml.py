# -*_ coding: utf-8 -*-
# __author__ = 'gaohao'
# last_modified: 20180408
from biocluster.core.exceptions import FileError
from biocluster.iofile import File
import os,re


class CgviewXmlFile(File):
    """
    只是对文件名称后缀做检查
    """
    def __init__(self):
        super(CgviewXmlFile, self).__init__()

    def check(self):
        super(CgviewXmlFile, self).check()
        if not os.path.exists(self.prop["path"]):
            raise FileError("%s文件不存在", variables=(self.prop["path"]), code="41400201")
        if not re.search(r'\.xml$',self.prop["path"]):
            raise FileError('文件必须以xml格式为后缀的文件！', code="41400202")



if __name__ == "__main__":
    a = CgviewXmlFile()
    a.check()

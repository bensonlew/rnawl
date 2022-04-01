# -*- coding: utf-8 -*-
# __author__ = 'wangbixuan'
from biocluster.iofile import File
from biocluster.core.exceptions import FileError

class GoListFile(File):
    """
    定义go_annotation输出的GO.list格式
    """
    def __init__(self):
        super(GoListFile,self).__init__()

    def get_info(self):
        super(GoListFile,self).get_info()
        with open(self.prop['path']) as f:
            go_info=f.read().split('\n')
            for record in go_info:
                if record!='':
                    go_record=record.split('\t')
                    for item in go_record[1].split(';'):
                        if not item.startswith('GO'):
                            raise FileError("存在错误的注释：%s", variables = (item), code = "45000601")

    def check(self):
        if super(GoListFile,self).check():
            self.get_info()
            return True
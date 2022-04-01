# -*_ coding: utf-8 -*-
# __author__ = 'zengjing'
# last_modified: 20201105
from biocluster.core.exceptions import FileError
from biocluster.iofile import File
import os


class SampleMergeListFile(File):
    """
    定义高通量数据拆分的输入sample_list文件的格式
    """
    def __init__(self):
        super(SampleMergeListFile, self).__init__()

    def check(self):
        super(SampleMergeListFile, self).check()
        if not os.path.exists(self.prop["path"]):
            raise FileError("%s文件不存在", variables=(self.prop["path"]), code="41800201")


if __name__ == "__main__":
    a = SampleMergeListFile()
    # a.set_path("")
    # a.check()
    # print a.prop["path"]

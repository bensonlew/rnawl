# -*_ coding: utf-8 -*-
# __author__ = 'HD'
# last_modified: 20181227
from biocluster.core.exceptions import FileError
from biocluster.iofile import File
import os


class FastqListFile(File):
    """
    检查第一列是样本名，第二列是样本fastq路径
    """
    def __init__(self):
        super(FastqListFile, self).__init__()

    def check(self):
        if super(FastqListFile, self).check():
            if not os.path.exists(self.prop["path"]):
                raise FileError("%s文件不存在" , variables=( self.prop["path"]), code="45500703")
            # else:
            #     self.check_len()
            return True

    # def check_len(self):
    #     with open(self.prop["path"], 'r') as r:
    #         if len(r.readline().split('\t')) != 2:
    #             raise FileError("%s文件列数不等于2！" % self.prop["path"])


if __name__ == "__main__":
    a = ListFileFile()
    a.set_path("/mnt/ilustre/users/sanger-dev/sg-users/zengjing/datasplit/ncRNA/nc_qc_after_list.txt")
    a.check()
    print a.prop["path"]

# -*_ coding: utf-8 -*-
# __author__ = 'zhaobinbin'
# last_modified: 20190117
from biocluster.core.exceptions import FileError
from biocluster.iofile import File
import os


class TagFileFile(File):

    def __init__(self):
        super(TagFileFile, self).__init__()

    def check(self):
        if super(TagFileFile, self).check():
            if not os.path.exists(self.prop["path"]):
                raise FileError("%s文件不存在" , variables=( self.prop["path"]), code="45500603")
            # else:
            #     self.check_len()
            return True

    # def check_len(self):
    #     with open(self.prop["path"], 'r') as r:
    #         if len(r.readline().split('\t')) != 2:
    #             raise FileError("%s文件列数不等于2！" % self.prop["path"])


if __name__ == "__main__":
    a = TagFileFile()
    a.set_path("/mnt/ilustre/users/sanger-dev/sg-users/zengjing/datasplit/ncRNA/nc_qc_after_list.txt")
    a.check()
    print a.prop["path"]
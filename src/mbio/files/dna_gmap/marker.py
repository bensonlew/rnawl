# -*- coding: utf-8 -*-

import os
from biocluster.iofile import File
from biocluster.core.exceptions import FileError


class MarkerFile(File):
    """
    定义marker格式文件
    hongdong@20180627
    属性：specimen_ids,marker文件里的所有样本list
          marker_num,所有的marker数目
    """

    def __init__(self):
        super(MarkerFile, self).__init__()

    def is_exists(self):
        if not os.path.isfile(self.path) or not os.path.exists(self.path):
            raise FileError("原始文件中不存在%s文件！", variables=(self.path), code="44800401")

    def check(self):
        if super(MarkerFile, self).check():
            self.is_exists()
            with open(self.path, "r") as f:
                lines = f.readlines()
                num = len(lines)
                if num < 2:
                    # raise FileError("文件:{}不得少于两行，请检查".format(self.path))
                    raise FileError("marker文件%s为空，请更改当前筛选条件，正常结束！", variables=(self.path), code="44800401")
                item = lines[0].strip().split("\t")
                specimen_ids = item[2:]
                # if len(item) <= 2:
                #     raise FileError("文件:{}不得少于两列，请检查！".format(self.path))
            self.set_property("specimen_ids", specimen_ids)
            self.set_property("marker_num", num-1)
            return True


if __name__ == '__main__':
    a = MarkerFile()
    a.set_path("/mnt/ilustre/users/sanger-dev/sg-users/zengjing/dna_gmap/cp_binmaker/Total.bin.marker")
    a.check()

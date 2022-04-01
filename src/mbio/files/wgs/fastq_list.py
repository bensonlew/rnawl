# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# modified 2018.04.03

import os
from biocluster.iofile import File
from biocluster.core.exceptions import FileError


class FastqListFile(File):
    """
    定义fastq.list文件
    第一列样本名称，第二列左端序列，第三列右端序列，第四列文库类型，列间以\t分隔
    """
    def __init__(self):
        super(FastqListFile, self).__init__()
        self.fastq_info = {}

    def check(self):
        if super(FastqListFile, self).check():
            self.get_info()
            return True

    def get_info(self):
        super(FastqListFile, self).get_info()
        self.get_fastq_info()
        self.set_property("fastq_info", self.fastq_info)

    def get_fastq_info(self):
        """
        得到字典self.fastq_info，样本及对应的双端序列
        """
        with open(self.path, "r") as f:
            for line in f:
                item = line.strip().split("\t")
                if len(item) < 3:
                    raise FileError("fastq.list里至少要有第一列样本，第二列左端序列，第三列右端序列", code="44500101")
                if item[0] not in self.fastq_info.keys():
                    self.fastq_info[item[0]] = {"l": item[1], "r": item[2]}
                else:
                    raise FileError("fastq.list里样本%s重复,请检查",variables=(item[0]), code="44500102")
        return self.fastq_info


if __name__ == "__main__":
    a = FastqListFile()
    a.set_path("/mnt/ilustre/users/sanger-dev/sg-users/zengjing/wgs/test_file/01.RawData/fastq.list")
    a.check()

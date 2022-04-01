# -*- coding: utf-8 -*-
# __author__ = "zengjing"
# last_modify: 20190124

import os
import re
from biocluster.iofile import File
from collections import defaultdict
from biocluster.core.exceptions import FileError


class DnaLibInfoFile(File):
    """
    定义DNA文库信息表
    表头：ProjectID\tLibID\tLibType\tSampleID\tEnzyme1\tEnzyme2\n
    """
    def __init__(self):
        super(DnaLibInfoFile, self).__init__()

    def check(self):
        super(DnaLibInfoFile, self).check()
        if not os.path.exists(self.prop["path"]):
            raise FileError("文件:{}不存在，请检查".format(self.prop["path"]))
        lib_exzyme, lib_type = defaultdict(list), defaultdict(list)
        lib_sample, pe_sample = [], []
        with open(self.prop["path"], "r") as f:
            lines = f.readlines()
            for line in lines[1:]:
                item = line.strip().split("\t")
                if len(item) < 4:
                    raise FileError("文件：{}每行长度必须大于等于4,请检查".format(self.prop["path"]))
                name = item[-1] + ":" + item[0] + ":" + item[1] + ":" + item[3]
                if name in lib_sample:
                    raise FileError("文库：{}中的样本：{}出现重复，请检查".format(name.split(":")[0], name.split(":")[1]))
                lib_sample.append(name)
                if re.search(r"PE", item[2]):
                    pe_sample.append(name)
                    lib_type[item[1]].append("PE")
                elif re.search(r"WGS", item[2]):
                    pe_sample.append(name)
                    lib_type[item[1]].append("WGS")
                elif re.search(r"RAD", item[2]):
                    if item[4] in lib_exzyme[item[1]]:
                        raise FileError("文库:{}中的酶: {}出现重复，请检查".format(item[1], item[4]))
                    lib_exzyme[item[1]].append(item[4])
                    lib_type[item[1]].append("RAD")
                elif re.search(r"GBS", item[2]):
                    enzyme = item[4] + ":" + item[5]
                    if enzyme in lib_exzyme[item[1]]:
                        raise FileError("文库:{}中的酶: {}-{}出现重复，请检查".format(item[1], item[4], item[5]))
                    lib_exzyme[item[1]].append(enzyme)
                    lib_type[item[1]].append("GBS")
                else:
                    raise FileError("文库:{}的类型不对，DNA中只能是PE/WGS/RAD/GBS".format(item[1]))
        for lib in lib_type.keys():
            if len(list(set(lib_type[lib]))) != 1:
                raise FileError("文库:{}的文库类型不唯一，请检查".format(lib))


if __name__ == "__main__":
    a = DnaLibInfoFile()
    a.set_path("/mnt/ilustre/users/sanger-dev/sg-users/zengjing/datasplit/datasplit_v2/dna_file/dna_lib_info.txt")
    a.check()

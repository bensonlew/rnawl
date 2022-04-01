# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# last_modify: 20180703

from biocluster.iofile import Directory
from biocluster.core.exceptions import FileError
import os


class TraitDirFile(Directory):
    """
    性状文件夹
    """
    def __init__(self):
        super(TraitDirFile, self).__init__()
        self.trait_file = {}

    def get_info(self):
        """
        获取文件属性
        """
        super(TraitDirFile, self).get_info()
        filelist = os.listdir(self.prop["path"])
        for file in filelist:
            if file.endswith(".txt"):
                pass
            else:
                raise FileError("文件夹中必须都是txt格式的文件", code="44801301")
            trait = file.split(".txt")[0]
            self.trait_file[trait] = os.path.join(self.prop["path"], file)
        self.set_property("trait_file", self.trait_file)

    def check(self):
        """
        检测文件是否满足要求，发生错误时应该触发FileError异常
        """
        if super(TraitDirFile, self).check():
            self.get_info()
            return True
        else:
            raise FileError("文件格式错误", code="44801302")


if __name__ == "__main__":
    a = TraitDirFile()
    a.set_path('/mnt/ilustre/users/sanger-dev/sg-users/zengjing/dna_gmap/qtl/new_file/trit_dir')
    a.check()

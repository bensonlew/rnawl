# -*- coding: utf-8 -*-
# __author__ = 'hongdong'
# modified 2018.06.27

import os
import re
from biocluster.iofile import Directory
from biocluster.core.exceptions import FileError


class MapCycleDirFile(Directory):
    """
    检查谱图评估的输入文件，如果是cp类型的，必须有*.sexAver.map， *.male.map， *.female.map，*.correct.loc,
    如果是nocp的时候*.out, *.correct.marker
    """
    def __init__(self):
        super(MapCycleDirFile, self).__init__()

    def dir_check(self):
        if 'path' in self.prop.keys() and os.path.isdir(self.prop['path']):
            n = 0
            for m in os.listdir(self.prop['path']):
                if not re.match(r".*\.pri\..*", m):
                    cp = re.match(r"(.*)\.sexAver\.map$", m)
                    if cp:
                        n += 1
                        self.is_file(os.path.join(self.prop['path'], "{}.male.map".format(cp.group(1))))
                        self.is_file(os.path.join(self.prop['path'], "{}.female.map".format(cp.group(1))))
                        self.is_file(os.path.join(self.prop['path'], "{}.correct.loc".format(cp.group(1))))
                    nocp = re.match(r'(.*)\.out$', m)
                    if nocp:
                        n += 1
                        self.is_file(os.path.join(self.prop['path'], "{}.correct.marker".format(nocp.group(1))))
            if n == 0:
                raise FileError("文件夹路径中没有图谱评估需要计算的文件!", code="44800301")
        else:
            raise FileError("文件夹路径不正确，请设置正确的文件夹路径--file error!", code="44800302")

    def is_file(self, file_path):
        """
        检查是否是文件是否存在
        :param file_path:
        :return:
        """
        if not os.path.isfile(file_path) or not os.path.exists(file_path):
            raise FileError("原始文件夹中不存在%s文件！", variables=(file_path), code="44800303")

    def check(self):
        if super(MapCycleDirFile, self).check():
            self.dir_check()
            return True

if __name__ == "__main__":
    a = MapCycleDirFile()
    a.set_path("/mnt/ilustre/users/sanger-dev/sg-users/xuanhongdong/gmap/05.map-cycle1-18-1-4")
    a.check()

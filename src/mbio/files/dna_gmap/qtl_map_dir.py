# -*- coding: utf-8 -*-
# __author__ = 'hongdong'
# modified 2018.06.27

import os
import re
from biocluster.iofile import Directory
from biocluster.core.exceptions import FileError


class QtlMapDirFile(Directory):
    """
    检查谱图评估的输入文件，如果是cp类型的，必须有*.sexAver.map， *.male.map， *.female.map，total.sexAver.loc,
    """
    def __init__(self):
        super(QtlMapDirFile, self).__init__()

    def dir_check(self):
        if 'path' in self.prop.keys() and os.path.isdir(self.prop['path']):
            cp_type = self.is_cp_type()
            nocp_type = self.is_nocp_type()
            if cp_type or nocp_type:
                pass
            else:
                raise FileError("文件夹路径: %s里的文件不正确，请设置正确的文件夹!" , variables=( self.prop["path"]), code="44801505")
        else:
            raise FileError("文件夹路径: %s路径不正确，请设置正确的路径!" , variables=( self.prop["path"]), code="44801506")

    def is_cp_type(self):
        """
        检查文件夹是否是F1类型的map文件夹
        :param file_path:
        :return:
        """
        male_map = os.path.join(self.prop['path'], "total.male.map")
        female_map = os.path.join(self.prop['path'], "total.female.map")
        sex_map = os.path.join(self.prop['path'], "total.sexAver.map")
        sex_loc = os.path.join(self.prop['path'], "total.sexAver.loc")
        cp_type = True
        if not os.path.exists(male_map):
            cp_type = False
        if not os.path.exists(female_map):
            cp_type = False
        if not os.path.exists(sex_map):
            cp_type = False
        if not os.path.exists(sex_loc):
            cp_type = False
        return cp_type

    def is_nocp_type(self):
        """
        检查文件夹是否是F2类型的map文件夹
        :return:
        """
        total_csv = os.path.join(self.prop['path'], "total.csv")
        total_map = os.path.join(self.prop['path'], "total.map")
        nocp_type = True
        if not os.path.exists(total_csv):
            nocp_type = False
        if not os.path.exists(total_map):
            nocp_type = False
        return nocp_type

    def check(self):
        if super(QtlMapDirFile, self).check():
            self.dir_check()
            return True

if __name__ == "__main__":
    a = QtlMapDirFile()
    a.set_path("/mnt/lustre/users/sanger/workspace/20181108/Single_nsanger_33014_QtlAnalysis_1108113127754418/QtlAnalysis/remote_input/total_map/evalutaion")
    a.check()

# -*- coding: utf-8 -*-
# __author__ = ysh
# last_modify：2019.04.15


import re, Bio, urllib2, regex, os
import subprocess
from biocluster.iofile import Directory
from collections import defaultdict
from biocluster.config import Config
from biocluster.core.exceptions import FileError


class MethyDirFile(Directory):
    """
    细菌基因组甲基化三代数据h5或bam文件夹
    """
    def __init__(self):
        super(MethyDirFile, self).__init__()

    def get_info(self):
        """
        获取文件属性
        :return:
        """
        if 'path' in self.prop.keys() and os.path.isdir(self.prop['path']):
            methy_sample_list,methy_sample_files,filedict,bam_dict= self.get_detail()
            self.set_property("methy_sample", methy_sample_list)
            self.set_property("methy_files", methy_sample_files)
            self.set_property("filedict", filedict)
            self.set_property("bam_dict", bam_dict)
        else:
            raise FileError("文件夹路径不正确，请设置正确的文件夹路径!")

    def check(self):
        """
        检测文件夹是否满足要求，不满足时触发FileError异常
        :return:
        """
        if super(MethyDirFile, self).check():
            self.get_info()
            return True

    def get_detail(self):
        """
        检测文件是否满足要求,发生错误时应该触发FileError异常
        head： sample	file	insert	length	size	lib
        :return:
        """
        filelist = os.listdir(self.prop['path'])
        methy_sample_list = []
        methy_sample_files = {}
        raw_list_path = ""
        if not len(filelist):
            raise FileError('甲基化文件夹为空，请检查确认')
        if "raw_list.txt" in filelist:
            raw_list_path = os.path.join(self.prop['path'], "raw_list.txt")
        elif  "list.txt" in filelist:
            raw_list_path = os.path.join(self.prop['path'], "list.txt")
        filedict = {}
        bam_dict = {}
        if raw_list_path:
            with open(raw_list_path, "rb") as inf:
                head = inf.next().strip().split("\t")
                #if head[0] != "sample" or head[1] != "file" or head[5] != "lib":
                #    raise FileError('list.txt文件格式有误，没有出现sample，file，lib字段')
                for line in inf:
                    line2 = line.strip().split("\t")
                    sample = line2[0]
                    files = line2[1]
                    lib = line2[5]
                    if lib == "pacbio" or lib =="Pacbio":
                        if "bax.h5" in files:
                            filedict[sample] = files
                            if not sample in methy_sample_list:
                                methy_sample_list.append(sample)
                                methy_sample_files[sample] = self.add_path(self.prop['path'], files)
                            else:
                                raise FileError('same sample {} has diff pacbio files'.format(sample))
                        elif ".bam" in files:
                            bam_dict[sample] = files
                            if not sample in methy_sample_list:
                                methy_sample_list.append(sample)
                                methy_sample_files[sample] = self.add_path(self.prop['path'], files)
                            else:
                                raise FileError('same sample {} has diff pacbio files'.format(sample))
        else:
            raise FileError('没有list文件，请输入！')
        return methy_sample_list,methy_sample_files,filedict,bam_dict

    def add_path(self, dir_path, files):
        files = files.split(",")
        new_list = []
        for each in files:
            new = os.path.join(dir_path, each)
            new_list.append(new)
        return ",".join(new_list)

    def get_suffix(self, files_str):
        '''
        检查三代数据后缀名是否为bax.h5或者bam
        '''
        filelist = files_str.split(",")
        sample_num =len(filelist)
        h5_num = 0
        bam_num = 0
        for filename in filelist:
            if "bax.h5" in filename:
                h5_num += 1
            elif "bam" in filename:
                bam_num += 1

        if h5_num > 0 :
            file_suffix = "h5"
        elif  bam_num  > 0:
            file_suffix = "bam"
        else:
            file_suffix = "fq"
                
        return file_suffix


if __name__ == "__main__":
    raw_dir = MethyDirFile()

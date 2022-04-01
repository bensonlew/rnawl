# -*- coding: utf-8 -*-
import os
import pandas as pd
from biocluster.iofile import File
from biocluster.core.exceptions import FileError
from biocluster.file import exists


class MethyFileFile(File):
    def __init__(self):
        super(MethyFileFile, self).__init__()

    def get_info(self):
        """
        获取文件属性
        :return:
        """
        if 'path' in self.prop.keys() and self.prop['path']:
            ref_samples,ref_dict= self.check_ref()
            self.set_property("ref_sample", ref_samples)
            self.set_property("ref_dict", ref_dict)
        else:
            raise FileError("文件夹路不正确，请设置正确的文件路径!")

    def check(self):
        """
        检测文件夹是否满足要求，不满足时触发FileError异常
        :return:
        """
        if super(MethyFileFile, self).check():
            self.get_info()
            return True

    def check_ref(self):
        '''
        检查ref_list file
        '''
        reflist = self.prop['path']
        ref_samples = []
        ref_dict = {}
        with open(reflist,"r") as inf:
            for line in inf:
                line = line.strip().split("\t")
                sample = line[0]
                ref_file = line[1]
                print sample
                if not ref_dict.has_key(sample):
                    ref_samples.append(sample)
                    ref_dict[sample] = ref_file
                else:
                    raise FileError('same sample {} has diff ref files'.format(sample))
        return ref_samples,ref_dict

    def check_ref_file(self):
        if super(MethyFileFile, self).check():
            if 'path' in self.prop.keys():
                reffile = self.prop['path']
                if not exists(reffile):
                    raise FileError('reference file {} has not exists'.format(reffile))
                else:
                    if "xml" in reffile:
                        file_suffix = "xml"
                    elif "fasta" in reffile:
                        file_suffix = "fasta"
                    else:
                        raise FileError("{} is not xml or fasta file".format(reffile))
                    self.set_property("ref_type", file_suffix)
        return True

    def get_ref_suffix(self, filelist):
        sample_num = len(filelist)
        xml_num = 0
        fasta_num = 0
        for filename in filelist:
            if "xml" in filename:
                file_suffix = "xml"
                xml_num += 1
            elif "fasta" in filename:
                file_suffix = "fasta"
                fasta_num += 1
            else:
                raise FileError("{} is not xml or fasta file".format(filename))
        ref_num = max(xml_num, fasta_num)
        if ref_num != sample_num:
           raise FileError("参考序列文件格式不一致")
        return file_suffix
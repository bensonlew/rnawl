# -*- coding: utf-8 -*-
# __author__ = gaohao
# last_modify：2018.03.16


import re, Bio, urllib2, regex, os
import subprocess
from biocluster.iofile import Directory
from collections import defaultdict
from biocluster.config import Config
from biocluster.core.exceptions import FileError
from mbio.files.sequence.fastq import FastqFile

'''
微生物基因组参数assemble文件夹检查
'''


class AsseDirFile(Directory):
    """
    定义raw_dir文件夹
    """
    def __init__(self):
        super(AsseDirFile, self).__init__()

    def get_info(self):
        """
        获取文件属性
        :return:
        """
        super(AsseDirFile, self).get_info()
        seqinfo = self.get_raw_info()
        self.set_property("sample_list", seqinfo[0])

    def check(self):
        """
        检测文件是否满足要求,发生错误时应该触发FileError异常
        :return:
        """
        filelist = os.listdir(self.prop['path'])
        asse_list_path = os.path.join(self.prop['path'], "list.txt")
        if not len(filelist):
            raise FileError('组装序列文件夹为空，请检查确认', code="41400101")
        if not os.path.exists(asse_list_path):
            raise FileError('组装序列文件夹的list.txt为不存在，请检查确认', code="41400102")
        files = os.listdir(self.prop['path'])
        with open(asse_list_path, "r") as f:
            lines = f.readlines()
            for line in lines:
                line2 = line.strip().split("\t")
                if len(line2) == 2 or len(line2) == 3:
                    pass
                else:
                    raise FileError('组装序列文件夹的list.txt不规范，请检查确认')
        for file in files:
            if re.search(r'.fasta',file) or re.search(r'.fna',file) or re.search(r'.fa',file):
                with open(self.prop['path'] + '/' + file,"r") as f, open(self.prop['path'] + '/' + file+'_new','w') as fw:
                    for lin in f:
                        lin = lin.rstrip('\n\r\t')
                        if not re.search('^>',lin):
                            if re.search('[^ATGCNatgcn]',lin):
                                #raise FileError('组装序列文件%s含有特殊字符！'%file)
                                lin = re.sub('[^atcgnATCGN]','N', lin)  #将兼并碱基替换成N
                        fw.write(lin+'\n')
                os.system('mv %s %s'%(self.prop['path'] + '/' + file+'_new', self.prop['path'] + '/' + file))

    def get_raw_info(self):
        """
        获取dir文件夹的样品信息
        :return: (sample_list)
        """
        sample_list ={}
        raw_list_path = os.path.join(self.prop['path'], "list.txt")
        with open(raw_list_path, "rb") as l:
            lines = l.readlines()
            for line in lines[1:]:
                line2 = line.strip().split("\t")
                sample_list[line2[0]] = line2[0]
        return sample_list

if __name__ == "__main__":
    raw_dir = AsseDirFile()
    #raw_dir.set_path("/mnt/lustre/users/sanger/workspace/20190314/Bacgenome_i-sanger_165214/remote_input/asse_dir/assemble")
    #test = raw_dir.get_info()
    #raw_dir.check()

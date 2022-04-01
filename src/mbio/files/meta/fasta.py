# -*- coding: utf-8 -*-
# __author__ = 'guoquan'
from biocluster.iofile import File
from collections import defaultdict
import re
import subprocess
from biocluster.config import Config
import os
from biocluster.core.exceptions import FileError
from Bio import SeqIO


class FastaFile(File):
    """
    定义Fasta文件， 需安装seqstat工具软件
    """

    def __init__(self):
        super(FastaFile, self).__init__()
        self.seqstat_path = os.path.join(Config().SOFTWARE_DIR, "bioinfo/seq/biosquid_1.9g+cvs20050121/bin/seqstat")
        # fasta与gff相关的属性
        self._seq_ids = []
        self._seq_obj = []
        self._seq_id_len_dic = {}
        self.seq_type = ''  # 有 prot nucl两种类型


    def get_info(self):
        """
        获取文件属性
        :return:
        """
        super(FastaFile, self).get_info()
        seqinfo = self.get_seq_info()
        self.set_property("file_format", seqinfo[0])
        self.set_property("seq_type", seqinfo[1])
        self.set_property("seq_number", seqinfo[2])
        self.set_property("bases", seqinfo[3])
        self.set_property("longest", seqinfo[4])
        self.set_property("shortest", seqinfo[5])


    def check(self):
        """
        检测文件是否满足要求,发生错误时应该触发FileError异常
        :return:
        """
        # print self.prop
        self.get_info()
        if super(FastaFile, self).check():
            if self.prop['file_format'] != 'FASTA':
                raise FileError("文件格式错误", code="42701301")
            if self.prop["seq_number"] < 1:
                raise FileError("应该至少含有一条序列", code="42701302")
        self.get_all_seq_name()
        return True


    def get_seq_info(self):
        """
        获取Fasta信息
        :return: (format,seq_type,seq_number,bases,longest,shortest)
        """
        try:
            subpro = subprocess.check_output(self.seqstat_path + " " + self.prop['path'], shell=True)
            result = subpro.split('\n')
            fformat = re.split(r':\s+', result[5])[1]
            seq_type = re.split(r':\s+', result[6])[1]
            seq_number = re.split(r':\s+', result[7])[1]
            bases = re.split(r':\s+', result[8])[1]
            shortest = re.split(r':\s+', result[9])[1]
            longest = re.split(r':\s+', result[10])[1]
            return fformat, seq_type, seq_number, bases, longest, shortest
        except subprocess.CalledProcessError:
            raise FileError("seqstat 运行出错！", code="42701303")


    def get_all_seq_name(self):
        seq_name = defaultdict(int)
        for seq in SeqIO.parse(self.prop["path"], "fasta"):
            seq_name[seq.id] += 1
        dup_list = list()
        for k in seq_name.iterkeys():
            if seq_name[k] > 1:
                dup_list.append(k)
            if len(dup_list) > 0:
                str_ = "; ".join(dup_list)
                raise FileError("序列名:%s在输入的fasta文件里面重复", variables=(str_), code="42701304")
        return seq_name


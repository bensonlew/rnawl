# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'
import re
import os
import shutil
import subprocess
from Bio import SeqIO
from biocluster.core.exceptions import FileError
from biocluster.iofile import Directory


class FastaDirFile(Directory):
    """
    定义fasta文件夹
    需要biopython
    """
    def __init__(self):
        """
        :param fastas: 不带路径的fastq的文件名集合
        :param fastas_full: 带路径的fastq的文件名集合
        """
        super(FastaDirFile, self).__init__()
        self.fastas = list()
        self.fastas_full = list()
        self.has_work_dir = False
        self.work_dir = ""

    def get_info(self):
        """
        获取文件夹属性
        """
        if 'path' in self.prop.keys() and os.path.isdir(self.prop['path']):
            self.set_property("fasta_number", self.get_fasta_number())
        else:
            raise FileError("文件夹路径不正确，请设置正确的文件夹路径!", code="44400501")
        self.get_check()

    def get_check(self):
        """
        检查
        :return:文件数目
        """
        list_txt = os.path.join(self.prop['path'], "list.txt")
        sample = []
        files = []
        if os.path.exists(list_txt):
            with open(list_txt, 'r') as f:
                lines = f.readlines()
                for i in lines[1:]:
                    lin = i.strip().split("\t")
                    if len(lin) != 2:
                        raise FileError('list.txt必须是两列！', code="44400512")
                    else:
                        files.append(lin[1])
                        sample.append(lin[0])
        filelist = os.listdir(os.path.join(self.prop['path']))
        if len(sample) != len(set(sample)):
            raise FileError('样品名称有一样的！', code="44400513")
        if len(files) != len(filelist)-1:
            raise FileError('list文件记录个数和文件夹下的个数不一致！', code="44400514")

    def get_full_info(self, work_path):
        """
        建立与这个fastq_dir相关的文件夹，获取全部的信息

        :param work_path: 工作文件夹的路径
        """
        self.make_work_dir(work_path)
        self.set_property("fasta_number", self.get_fasta_number())
        self.set_property("fasta_basename", self.fastas)

    def set_file_number(self, number):
        """
        设定文件中期望的fasta文件数，会与实际检测到的fasta做一个比较检验

        :param number: 设定的文件数
        """
        self.prop["expect_number"] = number

    def get_fasta_number(self):
        """
        获取文件夹下fasta的数目
        :return:文件数目
        """
        filelist = os.listdir(self.prop['path'])
        count = 0
        self.fastas = list()
        self.fastas_full = list()
        if not len(filelist):
            raise FileError('Fasta 序列文件夹为空，请检查确认', code="44400502")
        for file_ in filelist:
            if re.search(r'\.(fasta|fa|fna)$', file_):  # 增加fna后缀文件 by GHD @ 20180509
                count += 1
                self.fastas.append(file_)
                full_name = os.path.join(self.prop['path'], file_)
                self.fastas_full.append(full_name)
        return count


    def check(self):
        """
        检测文件夹是否满足要求，不满足时触发FileError异常
        :return:
        """
        if super(FastaDirFile, self).check():
            self.get_info()
            return True
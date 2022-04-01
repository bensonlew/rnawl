# -*- coding: utf-8 -*-
# __author__ = 'xuting'
import re
import os
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
            raise FileError("文件夹路径不正确，请设置正确的文件夹路径!", code="44000601")
        self.set_property("fasta_fullname", self.fastas_full)

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

    def make_work_dir(self, work_path):
        """
        创建临时文件夹

        :param work_path: 工作文件夹的路径
        """
        if not os.path.exists(work_path):
            os.mkdir(work_path)
        if os.path.isdir(work_path):  # 防止os.mkdir失败，做一次检测
            self.work_dir = work_path
            self.has_work_dir = True

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
            raise FileError('Fasta 序列文件夹为空，请检查确认', code="44000602")
        for file_ in filelist:
            if re.search(r'\.(fasta|fa|fna)$', file_):  # 增加fna后缀文件 by GHD @ 20180509
                count += 1
                self.fastas.append(file_)
                full_name = os.path.join(self.prop['path'], file_)
                self.fastas_full.append(full_name)
        return count

    def cat_fastas(self):
        """
        将所有的fasta文件合并到一起
        :return: 合并到一起的fasta文件路径
        """
        if not self.has_work_dir:
            raise FileError("还未建立工作路径！", code="44000603")
        cat_fasta = self.work_dir + "/cat_fasta.fasta"
        if os.path.exists(cat_fasta):
            os.remove(cat_fasta)
        os.mknod(cat_fasta)
        for fasta in self.fastas_full:
            try:
                cat_str = "cat " + fasta + " >> " + cat_fasta
                subprocess.check_call(cat_str, shell=True)
                return cat_fasta
            except subprocess.CalledProcessError:
                # error = "合并 " + fasta + u" 文件时出错"
                raise FileError("合并%s文件时出错", variables=(fasta), code="44000604")

    def cat_fastas_for_meta(self):
        """
        将所有的fasta按照一定的规则合并到一起
        规则: 检测fasta序列名里面有没有下划线_,有的话替换成单横杠-
              在序列名前加上文件名，与原序列名用下划线_隔开
              将所有的fasta合并到一起
        :return: 合并到一起的fasta文件的路径
        """
        if not self.has_work_dir:
            raise Exception("还未建立工作路径！", code="44000603")
        cat_fasta = self.work_dir + "/cat_meta.fasta"
        if os.path.exists(cat_fasta):
            os.remove(cat_fasta)
        os.mknod(cat_fasta)
        for fasta in self.fastas_full:
            basename = os.path.basename(fasta)
            sample_name = re.search(r'(.+)\.(fasta|fa)$', basename).group(1)
            with open(cat_fasta, "a") as f:
                for seq in SeqIO.parse(fasta, "fasta"):
                    new_id = str(sample_name) + '_' + str(seq.id)
                    f.write('>' + new_id + "\n")
                    f.write(str(seq.seq) + "\n")
        return cat_fasta

    def check(self):
        """
        检测文件夹是否满足要求，不满足时触发FileError异常
        :return:
        """
        if super(FastaDirFile, self).check():
            self.get_info()
            return True

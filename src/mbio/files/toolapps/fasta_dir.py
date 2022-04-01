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
        self.set_property("fasta_fullname", self.fastas_full)
        fasta_dir = self.check_unzip()
        self.set_property("fasta_dir", fasta_dir)
        self.get_check()

    def get_check(self):
        """
        检查
        :return:文件数目
        """
        list_txt = os.path.join(self.prop['path'], "fasta_dir", "list.txt")
        sample = []
        files = []
        if os.path.exists(list_txt):
            with open(list_txt, 'r') as f:
                lines = f.readlines()
                for i in lines:
                    lin = i.strip().split("\t")
                    if len(lin) != 2:
                        raise FileError('list.txt必须是两列！', code="44400512")
                    else:
                        files.append(lin[0])
                        sample.append(lin[1])
        filelist = os.listdir(os.path.join(self.prop['path'], "fasta_dir"))
        if len(sample) != len(set(sample)):
            raise FileError('样品名称有一样的！', code="44400513")
        if len(files) != len(filelist)-1:
            raise FileError('list文件记录个数和文件夹下的个数不一致！', code="44400514")
        else:
            new_filelist = []
            for file in filelist:
                file_path = os.path.join(self.prop['path'], "fasta_dir", file)
                if file not in new_filelist:
                    new_filelist.append(file_path)
            for i in files:
                if i not in new_filelist:
                    raise FileError('文件夹下{}不在list.txt文件中！'.format(i), code="44400515")
        self.prop['path'] = os.path.join(self.prop['path'], "fasta_dir")
        return self.prop['path']

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
            raise FileError('Fasta 序列文件夹为空，请检查确认', code="44400502")
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
            raise FileError("还未建立工作路径！", code="44400503")
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
                raise FileError("合并%s文件时出错", code="44400504")

    def cat_fastas_for_meta(self):
        """
        将所有的fasta按照一定的规则合并到一起
        规则: 检测fasta序列名里面有没有下划线_,有的话替换成单横杠-
              在序列名前加上文件名，与原序列名用下划线_隔开
              将所有的fasta合并到一起
        :return: 合并到一起的fasta文件的路径
        """
        if not self.has_work_dir:
            self.set_error("还未建立工作路径！", code="44400501")
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

    def check_unzip(self):
        """
        检验fasta文件夹是否是压缩格式，
        如果为压缩格式则进行解压，如果不是压缩格式则返回
        :return:
        """
        fasta_dir = os.path.join(self.prop['path'], "fasta_dir")
        old_list_file = os.path.join(self.prop['path'], "list.txt")
        list_file = os.path.join(self.prop['path'], "fasta_dir", 'list.txt')
        if os.path.exists(fasta_dir):
            shutil.rmtree(fasta_dir)
        filelist = os.listdir(self.prop['path'])
        os.mkdir(fasta_dir)
        if len(filelist) == 0:
            raise FileError('Fasta 序列文件夹为空，请检查确认', code="44400505")
        else:
            if os.path.exists(old_list_file):
                filelist.remove('list.txt')
                f = open(old_list_file, 'r')
                lines = f.readlines()
                sample_dict = {}
                for line in lines: ##获取原list文件的文件名称与样本名称的对应关系
                    line = line.strip().split("\t")
                    sample_name = line[1]
                    if sample_name not in sample_dict.keys():
                        sample_dict[line[0]] = sample_name
                with open(list_file, 'w') as w:
                    for file in filelist:
                        file_path = os.path.join(self.prop["path"], file)
                        if re.search(r'\.tar.gz$', file):
                            file_name = file.strip(".tar.gz")
                            unzip_cmd = "tar -zxvf " + file_path
                            try:
                                subprocess.check_output(unzip_cmd, shell=True)
                                old_path = os.path.join(self.prop["path"], file_name)
                                new_path = os.path.join(fasta_dir, file_name)
                                if os.path.exists(new_path):
                                    os.remove(new_path)
                                os.link(old_path, new_path)
                                if file in sample_dict.keys():
                                    w.write(new_path + "\t" + sample_dict[file] + "\n")
                            except subprocess.CalledProcessError:
                                raise FileError("%s:解压运行出错！", variables=(file), code="44400506")
                        elif re.search(r'\.gz$', file):
                            dir_path = fasta_dir
                            file_name = file.strip(".gz")
                            self.new_fasta = os.path.join(dir_path, file_name)
                            unzip_cmd = "zcat " +file_path + " > " + self.new_fasta
                            try:
                                subprocess.check_output(unzip_cmd, shell=True)
                                if file in sample_dict.keys():
                                    w.write(self.new_fasta + "\t" + sample_dict[file] + "\n")
                            except subprocess.CalledProcessError:
                                raise FileError("%s:解压运行出错！", variables=(file), code="44400507")
                        elif re.search(r'list', file):
                            pass
                        else:
                            self.new_fasta = os.path.join(fasta_dir, file)
                            os.link(file_path, self.new_fasta)
                            if file in sample_dict.keys():
                                w.write(self.new_fasta + "\t" + sample_dict[file] + "\n")
            else:
                with open(list_file, 'w') as w:
                    sample_list = []
                    for file in filelist:
                        file_path = os.path.join(self.prop["path"], file)
                        if re.search(r'\.tar.gz$', file):
                            file_name = file.strip(".tar.gz")
                            unzip_cmd = "tar -zxvf " + file_path
                            try:
                                subprocess.check_output(unzip_cmd, shell=True)
                                old_path = os.path.join(self.prop["path"], file_name)
                                new_path = os.path.join(fasta_dir, file_name)
                                if os.path.exists(new_path):
                                    os.remove(new_path)
                                os.link(old_path, new_path)
                                if re.search(r'\.fa$', file_name):
                                    sample_name = file_name.strip(".fa")
                                elif  re.search(r'\.fasta$', file_name):
                                    sample_name = file_name.strip(".fasta")
                                elif  re.search(r'\.fna$', file_name):
                                    sample_name = file_name.strip(".fna")
                                elif  re.search(r'\.ffn$', file_name):
                                    sample_name = file_name.strip(".ffn")
                                elif  re.search(r'\.fas$', file_name):
                                    sample_name = file_name.strip(".fas")
                                elif  re.search(r'\.faa$', file_name):
                                    sample_name = file_name.strip(".faa")
                                if sample_name not in sample_list:
                                    w.write(new_path + "\t" + sample_name + "\n")
                            except subprocess.CalledProcessError:
                                raise FileError("%s:解压运行出错！", variables=(file), code="44400508")
                        elif re.search(r'\.gz$', file):
                            dir_path = fasta_dir
                            file_name = file.strip(".gz")
                            self.new_fasta = os.path.join(dir_path, file_name)
                            unzip_cmd = "zcat " +file_path + " > " + self.new_fasta
                            try:
                                subprocess.check_output(unzip_cmd, shell=True)
                                if re.search(r'\.fa$', file_name):
                                    sample_name = file_name.strip(".fa")
                                elif  re.search(r'\.fasta$', file_name):
                                    sample_name = file_name.strip(".fasta")
                                elif  re.search(r'\.fna$', file_name):
                                    sample_name = file_name.strip(".fna")
                                else:
                                    raise FileError("%s:file没有正确的格式！", variables=(file), code="44400509")
                                if sample_name not in sample_list:
                                    w.write(self.new_fasta + "\t" + sample_name + "\n")
                            except subprocess.CalledProcessError:
                                raise FileError("%s:解压运行出错！", variables=(file), code="44400510")
                        else:
                            self.new_fasta = os.path.join(fasta_dir, file)
                            os.link(file_path, self.new_fasta)
                            if re.search(r"\.fa$", file):
                                sample_name = file.strip(".fa")
                            elif re.search(r'\.fasta$', file):
                                sample_name = file.strip(".fasta")
                            elif re.search(r'\.fna$', file):
                                sample_name = file.strip(".fna")
                            elif  re.search(r'\.ffn$', file_name):
                                    sample_name = file_name.strip(".ffn")
                            elif  re.search(r'\.fas$', file_name):
                                sample_name = file_name.strip(".fas")
                            elif  re.search(r'\.faa$', file_name):
                                sample_name = file_name.strip(".faa")
                            else:
                                raise FileError("%s:file没有正确的格式！", variables=(file), code="44400511")
                            if file not in sample_list:
                                w.write(self.new_fasta + "\t" + sample_name + "\n")
            return fasta_dir

    def check(self):
        """
        检测文件夹是否满足要求，不满足时触发FileError异常
        :return:
        """
        if super(FastaDirFile, self).check():
            self.get_info()
            return True

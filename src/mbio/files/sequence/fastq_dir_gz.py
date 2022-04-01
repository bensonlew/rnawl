# -*- coding: utf-8 -*-
# __author__ = 'xuting'
import re
import os
import subprocess
from biocluster.config import Config
from biocluster.iofile import Directory
from biocluster.core.exceptions import FileError
from mbio.files.sequence.fastq import FastqFile
from mbio.files.sequence.file_sample import FileSampleFile


class FastqDirFile(Directory):
    """
    定义fastq文件夹
    需要安装gzip
    需要安装fastoolkit，命令fastq_to_fasta用于将fastq转化为fasta
    """
    def __init__(self):
        """
        :param fastqs: 不带路径的fastq的文件名集合
        :param unzip_file: 带路径的fastq文件名的集合
        """
        super(FastqDirFile, self).__init__()
        self.fastq_to_fasta_path = os.path.join(Config().SOFTWARE_DIR, "bioinfo/seq/fastx_toolkit_0.0.14/fastq_to_fasta")
        self.fastqs = list()
        self.file_sample = dict()  # 文件名与样本名的对应，文件名为带绝对路径的全名
        self.is_convert = False
        self.has_unziped = False
        self.unzip_file = list()
        self.has_work_dir = False
        self.work_dir = ""
        self.has_list_file = False
        self.samples = list()
        self.se_repeat = False
        self.pe_repeat = False
        self.judge_the_gzs()

    def judge_the_gzs(self):
        list_txt = os.path.join(self.prop['path'], "list.txt")
        list_change = ''
        if os.path.exists(list_txt):
            with open(list_txt, 'r') as list_r:
                for line in list_r:
                    line = line.strip().split('\t')
                    if re.search(r'\.(fastq|fq)\.gz', line[0]):
                        ungz_name = re.search(r'(.+)\.(fastq|fq)\.gz$', line[0]).group(1)
                        new_fastq = os.path.join(self.work_dir, ungz_name + ".fastq")
                        try:
                            subprocess.check_call('gunzip -c ' + line[0] + " > " + new_fastq, shell=True)
                        except subprocess.CalledProcessError:
                            raise Exception("解压缩文件失败!")
                        list_change += new_fastq + '\t' + '\t'.join(line[1:]) + '\n'
                    else:
                        list_change += '\t'.join(line) + '\n'
            os.remove(list_txt)
            with open(list_txt, 'w') as list_w:
                list_w.write(list_change)

    def get_info(self):
        """
        获取文件夹属性
        """
        if 'path' in self.prop.keys() and os.path.isdir(self.prop['path']):
            self.set_property("fastq_number", self.get_fastq_number())
            self.set_property("fastq_basename", self.fastqs)
            self.set_property("has_list_file", self.has_list_file)
            self.set_property("samples", self.samples)
            self.set_property("file_sample", self.file_sample)
        else:
            raise FileError("文件夹路径不正确，请设置正确的文件夹路径!")

    def get_full_info(self, work_path):
        """
        建立与这个fastq_dir相关的文件夹，获取全部的信息(包括fastq解压)

        :param work_path: 工作文件夹的路径
        """
        self.get_info()
        self.make_work_dir(work_path)
        self.unzip_fastq()
        self.set_property("unzip_fastqs", list(set(self.unzip_file)))

    def get_fastq_number(self):
        """
        获取文件夹下fastq的数目
        :return:文件数目
        """
        list_txt = os.path.join(self.prop['path'], "list.txt")
        # test = os.path.join("/mnt/ilustre/users/sanger/test_xuting/otu/file_check/aaa.txt")
        # with open(test, "ab") as a:
        #     a.write(list_txt + "\n")
        if os.path.exists(list_txt):
            self.has_list_file = True
            filesample = FileSampleFile()
            filesample.set_path(list_txt)
            filesample.get_info()
            filesample.check()  # add 1 line by qiuping 20160722
            self.samples = filesample.prop["sample_names"]
            self.se_repeat = filesample.se_repeat
            self.pe_repeat = filesample.pe_repeat
            if not len(filesample.prop["file_names"]):
                raise FileError('Fastq 序列文件夹为空，请检查确认')
            for filename in filesample.prop["file_names"]:
                my_fastq = FastqFile()
                fq_path = os.path.join(self.prop['path'], filename)
                my_fastq.set_path(fq_path)
                my_fastq.get_info()
                sample_name = filesample.prop["file_sample"][filename]
                if my_fastq.check():
                    if filename not in self.fastqs:
                        self.fastqs.append(filename)
                        self.file_sample[fq_path] = sample_name
        else:
            filelist = os.listdir(self.prop['path'])
            if not len(filelist):
                raise FileError('Fastq 序列文件夹为空，请检查确认')
            for filename in filelist:
                if os.path.isdir(filename):
                    raise FileError('fastq文件夹中不应该存在文件夹：{}！！！'.format(filename))
                my_fastq = FastqFile()
                fq_path = os.path.join(self.prop['path'], filename)
                my_fastq.set_path(fq_path)
                my_fastq.get_info()
                sample_name = re.sub("\..+$", "", filename)
                if my_fastq.check():
                    self.fastqs.append(filename)
                    self.file_sample[fq_path] = sample_name

        return len(self.fastqs)

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

    def covert_to_fasta(self):
        """
        将所有的fastq转化为fasta文件，当fastq是gz格式的时候，先解压再转化
        :return:转化后的fasta文件夹地址
        """
        if not self.has_work_dir:
            raise Exception("还未建立工作路径！")
        if not self.is_convert:
            if not os.path.exists(os.path.join(self.work_dir, 'converted_fastas')):
                os.mkdir(os.path.join(self.work_dir, 'converted_fastas'))
            if self.has_unziped:
                for fastq in self.unzip_file:
                    fasta = re.search(r'(.+)\.(fastq|fq)', fastq).group(1)
                    fasta = os.path.join(self.work_dir, 'converted_fastas', os.path.basename(fasta) + ".fasta")
                    convert_str = (self.fastq_to_fasta_path + ' -Q 33' + ' -n -i '
                                   + fastq + ' -o ' + fasta)
                    mycmd = subprocess.Popen(convert_str, shell=True,
                                             stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
                    res = mycmd.communicate()
                    if mycmd.returncode == 0:
                        self.is_convert = True
                    else:
                        raise Exception('fastq转化fasta失败！\n' + convert_str + "\n" +
                                        res[0])
            else:
                raise Exception('文件还没有解压')
        return os.path.join(self.work_dir, 'converted_fastas')

    def unzip_fastq(self):
        """
        将压缩的fastq解压
        """
        if not self.has_work_dir:
            raise Exception("还未建立工作路径！")
        if not self.has_unziped:
            for fastq in self.prop["fastq_basename"]:
                fastq = os.path.join(self.prop['path'], fastq)
                if re.search(r'\.(fastq|fq)\.gz', fastq):
                    ungz_name = re.search(r'(.+)\.(fastq|fq)\.gz$', fastq).group(1)
                    new_fastq = os.path.join(self.work_dir, ungz_name + ".fastq")
                    try:
                        subprocess.check_call('gunzip -c ' + fastq + " > " + new_fastq, shell=True)
                        if new_fastq not in self.unzip_file:
                            self.unzip_file.append(new_fastq)
                    except subprocess.CalledProcessError:
                        raise Exception("解压缩文件失败!")
                else:
                    if fastq not in self.unzip_file:
                        new_fastq = os.path.join(self.work_dir, os.path.basename(fastq))  # added and edited by shijin
                        if not os.path.exists(new_fastq):
                            os.link(fastq, new_fastq)
                        self.unzip_file.append(new_fastq)
            self.has_unziped = True

    def check(self):
        """
        检测文件夹是否满足要求，不满足时触发FileError异常
        :return:
        """
        if super(FastqDirFile, self).check():
            self.get_info()
            return True

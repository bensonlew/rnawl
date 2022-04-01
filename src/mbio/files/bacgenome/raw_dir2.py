# -*- coding: utf-8 -*-
# __author__ = gaohao
# last_modify：2018.03.16


import re, regex, os
import subprocess
from biocluster.iofile import Directory
from biocluster.config import Config
from biocluster.core.exceptions import FileError
import pandas as pd
from biocluster.iofile import File
from mbio.files.sequence.fastq import FastqFile

'''
微生物基因组拼接流程输入参数rawdata文件检查
'''


class RawDir2File(Directory):
    """
    定义raw_dir文件夹
    """

    def __init__(self):
        super(RawDir2File, self).__init__()
        self.has_correct = False
        self.basic_info = {}  # 暂且不用
        self.samples = {}  # 验证后的样本，包含这个样本的拼接具体流程和基因组大小{samplename: (type, genomesize)}
        self.gz_split = {}  # 区分是否需要解压  {filepath: is zip or not }

    def get_info(self):
        """
        获取文件属性
        :return:
        """
        super(RawDir2File, self).get_info()
        seqinfo = self.get_raw_info()
        self.set_property("sample_list", seqinfo[0])
        # 要求可以判断出分析流程的类型: 纯二代/纯三代/混合

    def file_type_check(self, path, file_type):
        my_file = file_type()
        my_path = os.path.join(self.prop['path'], path)
        my_file.set_path(my_path)
        my_file.check()
        return my_file

    def check_sample_line(self, sample_data):
        # 获取只有一行数据的样品信息
        sample_task_type = None
        path_list = []
        sample_name = sample_data["Sample Name"].iloc[0]
        genome_size = 5 if sample_data["Genome Size(Mb)"].iloc[0] == "-" else sample_data["Genome Size(Mb)"].iloc[0]
        if sample_data.Library.iloc[0] == "PE":
            paths = sample_data["File Name"].iloc[0].split(",")
            for path in paths:
                my_fastq = self.file_type_check(path, FastqFile)
                if my_fastq.is_gz:
                    self.gz_split[path] = "gz"
                    # path_list.append({path: "gz"})
                else:
                    self.gz_split[path] = "not gz"
                    # path_list.append({path: "not gz"})
            sample_task_type = "draft"
            sample_info = {
                "task_type": sample_task_type,
                "insert_size": sample_data["Insert Size(bp)"].iloc[0],
                "read_length": sample_data["Reads Length(bp)"].iloc[0],
                "genome_size": genome_size,
                "library": sample_data["Library"].iloc[0],
                "path": path_list
            }
        elif sample_data.Library.iloc[0] == "MP":
            pass  # 不应该存在只有MP文库的数据
        else:
            paths = sample_data["File Name"].iloc[0].split(",")
            for path in paths:
                my_file = self.file_type_check(path, File)
                # path_list.append(path)
                if re.search('\.gz$', path) or re.search('\.gzip$', path):
                    self.gz_split[path] = "gz"
                elif re.search('\.tar\.gz$', path):
                    raise FileError("不支持tar.gz格式文件")
                elif re.search('\.bas.h5', path):
                    raise FileError("流程中的h5文件中不需要bas.h5文件: %s" % path)
                else:
                    self.gz_split[path] = "not gz"
            sample_task_type = "complete"
            sample_info = {
                "task_type": sample_task_type,
                "insert_size": "-",
                "read_length": "-",
                "genome_size": genome_size,
                "library": sample_data["Library"].iloc[0],
                "path": path_list
            }
        if sample_task_type:
            self.basic_info[sample_name] = {
                "library_name": sample_data["Library Name"].iloc[0],
                "sample_info": sample_info
            }
            self.samples[sample_name] = (sample_task_type, genome_size)
        return sample_task_type

    def check_sample_data(self, sample, sample_data):
        # 获取有多行数据的样品信息
        self.basic_info[sample] = {}
        sample_library_set = set(sample_data.Library)
        genome_size_set = set(sample_data["Genome Size(Mb)"]) - set(["-"])
        if len(genome_size_set) > 1:
            error_str = ",".join(list(genome_size_set))
            raise FileError("样品%s有多个Genome Size: %s" % (sample, error_str))
        elif len(genome_size_set) == 1:
            genome_size = genome_size_set.pop()
        else:
            genome_size = 5
        if set(["Pacbio", "Nanopore"]) - sample_library_set == set(["Pacbio", "Nanopore"]):
            for file_name in sample_data["File Name"]:
                for path in file_name.split(","):
                    my_fastq = self.file_type_check(path, FastqFile)
                    if my_fastq.is_gz:
                        # path_list.append({path: "gz"})
                        self.gz_split[path] = "gz"
                    else:
                        # path_list.append({path: "not gz"})
                        self.gz_split[path] = "not gz"
            sample_task_type = "draft"
        else:
            for file_name in sample_data["File Name"]:
                for path in file_name.split(","):
                    my_file = self.file_type_check(path, File)
                    if re.search('\.gz$', path) or re.search('\.gzip$', path):
                        self.gz_split[path] = "gz"
                    elif re.search('\.tar\.gz$', path):
                        raise FileError("不支持tar.gz格式文件")
                    else:
                        self.gz_split[path] = "not gz"
            sample_task_type = "complete"
            sample_info = {
                "task_type": sample_task_type,
                "insert_size": "-",
                "read_length": "-",
                "genome_size": genome_size,
            }
        self.samples[sample] = (sample_task_type, genome_size)
        self.basic_info[sample]["task_type"] = sample_task_type
        return sample_task_type

    def check(self):
        """
        检测文件是否满足要求,发生错误时应该触发FileError异常
        :return:
        """
        filelist = os.listdir(self.prop['path'])
        raw_list_path = os.path.join(self.prop['path'], "list.txt")
        # sample_list = {}
        if not len(filelist):
            raise FileError('原始序列文件夹为空，请检查确认', code="41400801")
        if not os.path.exists(raw_list_path):
            raise FileError('原始序列文件夹的list.txt为不存在，请检查确认', code="41400802")
        data = pd.read_table(raw_list_path).fillna("-")
        # if set(data.columns) !=  set(['Library Name', 'File Name', 'Library', 'Genome Size(Mb)', 'Insert Size(bp)', 'Reads Length(bp)', 'Sample Name']):
        if set(data.columns) != set(['lib_name', 'file', 'lib', 'size', 'insert', 'length', 'sample']):
            if set(data.columns) !=  set(['Library Name', 'File Name', 'Library', 'Genome Size(Mb)', 'Insert Size(bp)', 'Reads Length(bp)', 'Sample Name']):
                raise FileError('原始序列文件夹的list.txt表头不符合规范，请检查确认')
        else:
            data.columns = ['Sample Name', 'File Name', 'Insert Size(bp)', 'Reads Length(bp)', 'Genome Size(Mb)', 'Library', 'Library Name']
        library_set = set(data.Library)
        task_type_confirm = False
        samples = set(data["Sample Name"])
        if not library_set <= set(['PE', 'MP', 'Pacbio', 'Nanopore']):
            error_str = ",".join(list(library_set - set(['PE', 'MP', 'Pacbio', 'Nanopore'])))
            raise FileError('Library必须为PE/MP/Pacbio/Nanopore中的一个，不能为：%s' % error_str)
        elif set(["Pacbio", "Nanopore"]) - library_set == set(["Pacbio", "Nanopore"]):
            self.set_property("task_type", "draft")
            task_type_confirm = True
        elif set(["MP", "PE"]) - library_set == set(["MP", "PE"]):
            self.set_property("task_type", "chr")
            task_type_confirm = True
        for sample in samples:
            sample_data = data[data["Sample Name"] == sample]
            if len(sample_data) > 1:
                if self.check_sample_data(sample, sample_data) == "draft" and not task_type_confirm:
                    '''
                    此处判断时task_type可能为"chr"或者"mix",两者肯定是有三代
                    数据的，但对应于一个样品来说，可能不含三代数据，所以如果这个样品是纯二代拼接的，
                    则可以断定这个总流程类型是"mix"
                    '''
                    self.set_property("task_type", "mix")
                    task_type_confirm = True
            else:
                if self.check_sample_line(sample_data) == "draft" and not task_type_confirm:
                    '''
                    道理同上
                    '''
                    self.set_property("task_type", "mix")
                    task_type_confirm = True
        if not task_type_confirm:  # 如果所有样品查看后仍没有发现扫描图拼接，则可以断定这个任务就是完成图拼接
            self.set_property("task_type", "chr")  # draft为纯扫描图，chr为纯完成图，mix为兼有

    def get_raw_info(self):
        """
        获取dir文件夹的样品信息
        :return: (sample_list)
        """
        sample_list = {}
        raw_list_path = os.path.join(self.prop['path'], "list.txt")
        # with open(raw_list_path, "rb") as l:
        #     lines = l.readlines()
        #     for line in lines[1:]:
        #         line2 = line.strip().split("\t")
        #         if len(line2) == 6:
        #             sample_list[line2[0]] = line2[0]
        sample_list = ["1"]
        return sample_list

    def list_correcct(self):
        """
        验证list文件，对于文库名称误写的情况进行校正，不能校正的样品排除掉
        校正成功的加入到样品的详细属性中
        :return:
        """
        self.has_correct = True
        pass

    def read_sample_info(self, path):
        """
        从已验证的样品信息中获取,便于工作流中通过修改中间样本属性文件调试流程
        :param path:
        :return:
        """
        self.has_correct = True

    def write_sample_info(self, outpath):
        pass

if __name__ == "__main__":
    raw_dir = RawDir2File()

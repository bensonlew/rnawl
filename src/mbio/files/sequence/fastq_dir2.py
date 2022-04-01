# -*- coding: utf-8 -*-
import os
from collections import defaultdict
from biocluster.iofile import Directory
from biocluster.core.exceptions import FileError
from mbio.files.sequence.fastq import FastqFile
from mbio.files.sequence.fastq_dir import FastqDirFile
from mbio.files.sequence.file_sample import FileSampleFile


class FastqDir2File(FastqDirFile):
    """
    定义fastq文件夹
    需要安装gzip
    """
    def __init__(self):
        """
        :param fastqs: 不带路径的fastq的文件名集合
        :param unzip_file: 带路径的fastq文件名的集合
        """
        super(FastqDir2File, self).__init__()
        self.fqs = []
        self.gzipped = ''
        self.sample_info = defaultdict(dict)  # 样本及对应的方向，文件路径
        self.sups = 0  # 存在补测的样本个数

    def fastq_check(self):
        for fq in self.fqs:
            try:
                fq.check()
            except FileError as e:
                raise FileError("{} fastq格式检查失败".format(os.path.basename(fq.path)))
            except Exception as e:
                raise Exception("{} 检查失败 {}".format(os.path.basename(fq.path)))

    def list_check(self):
        list_file = os.path.join(self.path, "list.txt")
        if not os.path.exists(list_file):
            raise FileError("list.txt文件不存在")
        try:
            list_check = FileSampleFile()
            list_check.set_path(list_file)
            list_check.check()
        except FileError as e:
            raise FileError("{} fastq格式检查失败 {}".format(list_file, e))
        list_check.check()
        file_list = list_check.get_list()
        sup_seq = {}
        not_gzip = []
        for sample, info in file_list.items():
            if isinstance(info, str):
                files_path = [os.path.join(self.path, f) for f in info.split(',')]
                self.sample_info[sample]['s'] = files_path
                for one in files_path:
                    fastq = FastqFile()
                    fastq.set_path(one)
                    self.fqs.append(fastq)
                    if not fastq.is_gzip:
                        not_gzip.append(os.path.basename(one))
            else:
                for di in info:
                    files_path = [os.path.join(self.path, f) for f in info[di].split(',')]
                    if len(files_path) > 1:
                        sup_seq[sample] = True
                    self.sample_info[sample][di] = sorted(files_path)
                    for one in files_path:
                        fastq = FastqFile()
                        fastq.set_path(one)
                        self.fqs.append(fastq)
                        if not fastq.is_gzip:
                            not_gzip.append(os.path.basename(one))
        if len(not_gzip) > 0:
            raise FileError("下列文件gzip压缩格式不正确:\n{}".format(not_gzip))
        self.sups = len(sup_seq)

    def check(self):
        """
        检测文件夹是否满足要求，不满足时触发FileError异常
        :return:
        """
        super(FastqDir2File, self).check()
        self.get_info()
        self.list_check()
        return True

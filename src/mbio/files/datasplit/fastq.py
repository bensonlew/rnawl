# -*- coding: utf-8 -*-
# __author__ = 'xuting'
import os
import re
import gzip
import subprocess
from biocluster.iofile import File
from biocluster.core.exceptions import FileError
from biocluster.config import Config


class FastqFile(File):
    """
    定义Fastq文件，
    需要安装gzip
    需要安装fastoolkit，命令fastq_to_fasta用于将fastq转化为fasta
    需要安装seqstat软件，版本1.9g(Oct 2002),软件用法:  seqstat <fasta文件>，根据stdout的输出，统计fasta的信息
    :param _fastaname:由fastq文件转化而来的fasta文件
    :param _filename: 若fastq为gz格式,为解压后fastq文件，若fastq不是gz格式，等于prop["path"]
    :param is_convert: 是否已经转化成fasta
    """
    def __init__(self):
        super(FastqFile, self).__init__()
        #self.seqstat_path = os.path.join(Config().SOFTWARE_DIR, "seqs/seqstat")
        self.seqstat_path = os.path.join(Config().SOFTWARE_DIR, "bioinfo/seq/biosquid_1.9g+cvs20050121/bin/seqstat") #修改软件路径 konghualei 20170124
        self._fastaname = ""
        self._filename = ""
        #self.fastq_to_fasta_path = os.path.join(Config().SOFTWARE_DIR, "fastxtoolkit/bin/fastq_to_fasta")
        self.fastq_to_fasta_path = os.path.join(Config().SOFTWARE_DIR, "bioinfo/seq/fastx_toolkit_0.0.14/fastq_to_fasta") # konghualei 20170124
        self.is_convert = False
        self.has_sample_info = False
        self.samples = list()
        self.is_gz = False

    @property
    def is_gzip(self):
        """
        依据文件后缀名检测是是gz个是压缩文件
        """
        if re.search('\.tar\.gz$', self.prop['path']):
            raise FileError("不支持tar.gz格式文件")
        if re.search('\.gz$', self.prop['path']) or re.search('\.gzip$', self.prop['path']):
            return True
        else:
            return False

    def get_info(self):
        """
        获取文件的基本属性
        """
        super(FastqFile, self).get_info()
        format_ = self.check_format()
        self.set_property("is_gz", self.is_gz)
        self.set_property("has_sample_info", self.has_sample_info)
        self.set_property("is_fqformat", format_)

    def get_full_info(self, work_path):
        """
        获取文件属性
        当fastq是gz格式的时候,解压到work_path
        将fastq转化成fasta，并统计相关信息
        :param work_path:工作文件夹路径
        """
        self.get_info()
        self._prepare(work_path)
        self.convert_to_fasta()
        seqinfo = self.get_seq_info()
        self.set_property("unzip", self.unzipfile)
        self.set_property("fasta", self.fastaname)
        self.set_property("fasta_fomat", seqinfo[0])
        self.set_property("seq_number", seqinfo[1])
        self.set_property("bases", seqinfo[2])
        self.set_property("longest", seqinfo[3])
        self.set_property("shortest", seqinfo[4])

    def check(self):
        """
        检测文件是否满足要求,发生错误时应该触发FileError异常
        """
        if super(FastqFile, self).check():
            self.get_info()
            return self.prop["is_fqformat"]

    def check_format(self):
        """
        检测文件是否满足要求,发生错误时应该触发FileError异常
        :return: bool
        """
        self.is_gz = self.is_gzip
        if self.is_gz:
            with gzip.open(self.prop['path'], 'rb') as f:
                line1 = f.next()
                line = f.next()
                line = f.next()
                line = f.next()
                try:
                    line5 = f.next()
                    if not (re.search(r'^@', line1) and re.search(r'^@', line5)):
                        raise FileError("非压缩后的fastq格式文件")
                    myline1 = re.split('_', line1)
                    myline2 = re.split('_', line5)
                    if len(myline1) > 1 and len(myline2) > 1:
                        self.has_sample_info = True
                except:
                    print "{}只有一条序列".format(self.prop['path'])
        else:
            with open(self.prop['path'], 'r') as r:
                line = r.next()
                if not re.search(r'^@', line):
                    raise FileError("fastq文件格式错误")
                myline1 = re.split('_', line)
                try:
                    line = r.next()
                    line = r.next()
                    line = r.next()
                    line = r.next()
                    if not re.search(r'^@', line):
                        raise FileError("fastq文件:{}格式错误".format(self.path))
                    myline2 = re.split('_', line)
                    if len(myline1) > 1 and len(myline2) > 1:
                        self.has_sample_info = True
                except:
                    print "{}只有一条序列".format(self.prop['path'])
        return True

    def check_content(self):
        """
        遍历一个fastq文件，查看是否符合fastq的规范， 并获得所有的样本名
        可能需要较长的时间
        """
        if self.is_gz:
            with gzip.open(self.prop['path'], 'rb') as f:
                count = 2
                for line in f:
                    count = count + 4
                    line = line.rstrip("\r\n")
                    if not re.search(r'^@', line):
                        raise Exception("未检测到@，非压缩后的fastq格式文件")
                    line = re.sub("^@", "", line)
                    line = re.split('\s+', line)[0]
                    line = re.split('_', line)
                    line.pop(-1)
                    sp_name = "_".join(line)
                    if sp_name not in self.samples:
                        self.samples.append(sp_name)
                    line2 = f.next()
                    f.next()
                    line4 = f.next()
                    if len(line2) != len(line4):
                        raise Exception("第{}行碱基的长度与它的质量文件的长度不相等".format(str(count)))
            self.set_property("samples", self.samples)
        else:
            with open(self.prop['path'], 'rb') as f:
                count = 2
                for line in f:
                    line = line.rstrip("\r\n")
                    if not re.search(r'^@', line):
                        raise Exception("未检测到@，非fastq格式文件 第{}行缺少@".format(count))
                    line = re.sub("^@", "", line)
                    line = re.split('\s+', line)[0]
                    line = re.split('_', line)
                    line.pop(-1)
                    sp_name = "_".join(line)
                    if sp_name not in self.samples:
                        self.samples.append(sp_name)
                    try:
                        line2 = f.next().rstrip()
                        f.next()
                        line4 = f.next().rstrip()
                    except Exception:
                        raise Exception("fastq文件不完整，请查看fastq文件的最后")
                    if len(line2) != len(line4):
                        raise Exception("第{}行碱基的长度与它的质量文件的长度不相等, 碱基长度为{}, 质量长度为{}".format(str(count), str(len(line2)), str(len(line4))))
                    count = count + 4
            self.set_property("samples", self.samples)
        return self.samples

    def _prepare(self, work_path):
        """
        为获取序列的信息做准备
        生成临时文件夹，当输入的文件是gz格式时，解压到tmp里
        """
        self.unzipfile = self.prop['path']
        basename = os.path.basename(self.prop['path'])
        self.fastaname = os.path.join(work_path, basename + ".fasta")
        if self.is_gz:
            basename = re.search(r'(.+)\.gz', basename).group(1)
            self.unzipfile = os.path.join(work_path, basename)
            try:
                subprocess.check_call('gunzip -c ' + self.prop['path'] + "> " + self.unzipfile)
            except subprocess.CalledProcessError:
                raise FileError("非标准格式的gz文件！")

    def convert_to_fasta(self):
        """
        将fastq转化成fasta
        :return: 转化后的fasta的文件
        """
        if not self.is_convert:
            try:
                convert_str = (self.fastq_to_fasta_path + ' -n -i '
                               + self.unzipfile + ' -o ' + self.fastaname)
                subprocess.check_call(convert_str, shell=True)
                self.is_convert = True
            except subprocess.CalledProcessError:
                raise Exception('fastq转化fasta失败！')
        return self.fastaname

    def get_seq_info(self):
        """
        获取Fasta信息
        :return: (format,seq_number,bases,longest,shortest)
        """
        try:
            subpro = subprocess.check_output(self.seqstat_path + " " + self.fastaname, shell=True)
            result = subpro.split('\n')
            seq_type = re.split(r':\s+', result[6])[1]
            seq_number = re.split(r':\s+', result[7])[1]
            bases = re.split(r':\s+', result[8])[1]
            shortest = re.split(r':\s+', result[9])[1]
            longest = re.split(r':\s+', result[10])[1]
            return seq_type, seq_number, bases, shortest, longest
        except subprocess.CalledProcessError:
            raise Exception("seqstat 运行出错！")

if __name__=="__main__":
    data=FastqFile()
    data.set_path("/mnt/ilustre/users/sanger-dev/sg-users/konghualei/ref_rna/tools/kallisto/fq_dir/fastq/trim1.fq")
    data._prepare("/mnt/ilustre/users/sanger-dev/sg-users/konghualei/ref_rna/tools/kallisto/fq_dir/fastq/")
    data.convert_to_fasta()
    d1=data.get_seq_info()
    print d1
    #a1=data.check_format()
    #data.get_info()

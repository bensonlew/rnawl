# -*- coding: utf-8 -*-
# __author__ = ''
# __modify__ = '20201216'


import os
import re
import pandas as pd
import time
import shutil
from collections import defaultdict
from biocluster.module import Module
from biocluster.core.exceptions import OptionError
from mbio.files.sequence.file_sample import FileSampleFile
from mbio.packages.metagbin.common_function import link_dir


class BacgenomeFastpModule(Module):
    """
    fastp的质控module
    """

    def __init__(self, work_id):
        super(BacgenomeFastpModule, self).__init__(work_id)
        option = [
            {"name": "sample_path", "type": "infile","format": "sequence.fastq_dir"},
            {'name': "sample_info", "type": "string"},
            {"name": "clean_list", "type": "outfile", "format": "bacgenome.list_file"},
            {"name": "fq_type", "type": "string", "default": "PE"},  # PE OR SE
            {"name": "qualified_quality_phred", "type": "string", "default": "20"},
            # -q,一个碱基合格的质量值,默认表示phred质量> = Q是合格的。
            {"name": "length_required", "type": "string", "default": "50"},  # -l,长度过滤参数，比此值短的读取将被丢弃
            {"name": "cut_by_quality5", "type": "string", "default": "20"},  # -5,根据前面(5 ")的质量，允许每个读切割，默认是禁用的
            {"name": "cut_by_quality3", "type": "string", "default": "20"},  # -3,根据后面(3 ")的质量，允许每个读切割，默认是禁用的
            {"name": "cut_mean_quality", "type": "string", "default": "20"},  # -M,在滑动窗口的基础上的平均质量低于切割质量将被切割，默认是Q20
            {"name": "n_base_limit", "type": "string", "default": "0"},  # -n,如果reads的碱基数大于该值，那么这个reads就被丢弃了
            {"name": "compression", "type": "string", "default": "5"},  # -z,gzip输出的压缩级别(1 ~ 9). 1是最快的，9是最小的
            {"name": "thread", "type": "string", "default": "8"},  # -w,线程数
        ]
        self.add_option(option)
        self.sample_path = defaultdict(list)
        self.sample_info = {}
        self.read_info = {}
        self.qc_tools = []

    def check_options(self):
        """
        检查参数
        :return:
        """
        if not self.option('sample_path').is_set:
            raise OptionError('必须输入样本文件夹对应的路径信息', code="21400601")
        row_num = len(open(os.path.join(self.option("sample_path").prop['path'],"list.txt"), "r").readline().split())
        if row_num != 3:
            raise OptionError("PE序列list文件应该包括文件名、样本名和左右端说明三列", code="21400602")
        if not self.option('sample_info'):
            raise OptionError('必须输入样本信息文件，包括样本，文库类型、插入片段长度三列', code="21400603")
        return True

    def run(self):
        super(BacgenomeFastpModule, self).run()
        self.get_info()
        time.sleep(2)
        self.run_qc()

    def run_qc(self):
        """
        将文件夹中的文件按样本名称去质控
        :return:
        """
        self.samples = self.get_list()
        self.logger.info("{}".format(self.samples))
        fastq_path = self.option("sample_path").prop['path']
        for key in self.samples.keys():
            fq_l = os.path.join(fastq_path, self.samples[key]["l"])
            fq_r = os.path.join(fastq_path, self.samples[key]["r"])
            self.qc = self.add_tool("bacgenome.fastp")
            opts = ({
                "fq1": fq_l,
                "fq2": fq_r,
                "qualified_quality_phred": self.option("qualified_quality_phred"),
                "length_required": self.option("length_required"),
                "cut_mean_quality": self.option("cut_mean_quality"),
                "n_base_limit": self.option("n_base_limit"),
                "thread": self.option("thread"),
                "cut_by_quality3": self.option("cut_by_quality3"),
                "cut_by_quality5": self.option("cut_by_quality5"),
                "sample_name": key
            })
            self.qc.set_options(opts)
            self.qc_tools.append(self.qc)
        if len(self.qc_tools)> 1:
            self.on_rely(self.qc_tools, self.set_output)
        else:
            self.qc_tools[0].on('end', self.set_output)
        if self.qc_tools:
            for tool in self.qc_tools:
                tool.run()
        else:
            self.set_error("tool 队列是空的，无法进行后续的分析，流程终止")

    def get_list(self):
        list_path = os.path.join(self.option("sample_path").prop["path"], "list.txt")
        file_sample = FileSampleFile()
        file_sample.set_path(list_path)
        samples = file_sample.get_list()
        return samples

    def set_output(self):
        """
        将结果文件连接到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        fastq_dir = os.path.join(self.output_dir)
        with open(os.path.join(fastq_dir, "list.txt"), "w") as w:
            for obj in self.qc_tools:
                for file in os.listdir(obj.work_dir):
                    if re.search(r'\.json', file):
                        s = str(file.split('.')[0])
                        fq1_path = obj.work_dir + "/" + s + ".clean.1.fastq"
                        fq2_path = obj.work_dir + "/" + s + ".clean.2.fastq"
                        unpaired_path = obj.work_dir + "/" + s + ".unpaired.fastq"
                        fq1_name = s + ".clean.1.fastq"
                        fq2_name = s + ".clean.2.fastq"
                        unpaired_name = s + ".unpaired.fastq"
                        fq1_path_new = os.path.join(fastq_dir, fq1_name)
                        fq2_path_new = os.path.join(fastq_dir, fq2_name)
                        unpaired_path_new = os.path.join(fastq_dir, unpaired_name)
                        if os.path.exists(fq1_path_new):
                            os.remove(fq1_path_new)
                        if os.path.exists(fq2_path_new):
                            os.remove(fq2_path_new)
                        if os.path.exists(unpaired_path_new):
                            os.remove(unpaired_path_new)
                        os.link(fq1_path, fq1_path_new)
                        os.link(fq2_path, fq2_path_new)
                        os.link(unpaired_path, unpaired_path_new)
                        w.write(fq1_name + "\t" + s + "\t" + "l" + "\n")
                        w.write(fq2_name + "\t" + s + "\t" + "r" + "\n")
                        w.write(unpaired_name + "\t" + s + "\t" + "s" + "\n")
        self.option('clean_list').set_path(self.output_dir + "/list.txt")
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [",", "", "结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            ["", "", ""]
        ])
        super(BacgenomeFastpModule, self).end()

    def linkdir(self, dirpath, dirname):
        """
        link一个文件夹下的所有文件到本module的output目录
        :param dirpath: 传入文件夹路径
        :param dirname: 新的文件夹名称
        :return:
        """
        newdir = os.path.join(self.output_dir, dirname)
        if not os.path.exists(newdir):
            os.mkdir(newdir)
        allfiles = os.listdir(dirpath)
        oldfiles = [os.path.join(dirpath, i) for i in allfiles]
        newfiles = [os.path.join(newdir, i) for i in allfiles]
        for newfile in newfiles:
            if os.path.exists(newfile):
                if os.path.isfile(newfile):
                    os.remove(newfile)
                else:
                    os.system('rm -r %s' % newfile)
        for i in range(len(allfiles)):
            if os.path.isfile(oldfiles[i]):
                os.link(oldfiles[i], newfiles[i])
            elif os.path.isdir(oldfiles[i]):
                os.link(oldfiles[i], newdir)

    def get_info(self):
        """
        获得样本对应的路径信息
        :return:
        """
        with open(self.option("sample_info"))as f:
            lines = f.readlines()
            for line in lines[0:]:
                tmp = line.strip().split('\t')
                self.sample_info[tmp[0]] = tmp[1]  # 样本的插入片段长度
                self.read_info[tmp[0]] =tmp[2]  #read长度
        with open(os.path.join(self.option("sample_path").prop['path'],"list.txt"))as fr:
            for line in fr:
                # self.logger.info(self.sample_path)
                tmp = line.strip().split('\t')
                if tmp[1] in self.sample_info.keys():
                    if tmp[1] in self.sample_path.keys():
                        if tmp[2] == 'l':
                            self.sample_path[tmp[1]].insert(0, tmp[0])
                        else:
                            self.sample_path[tmp[1]].append(tmp[0])
                    else:
                        self.sample_path[tmp[1]].append(tmp[0])
                else:
                    self.set_error('需要质控的序列样本%s没有相关的样本信息，请核实！', variables=(tmp[1]), code="21400601")
            for key in self.sample_path.keys():
                if len(self.sample_path[key]) > 2:
                    self.set_error('需要质控的序列样本%s有重名，请改样本名或分开质控！', variables=(key), code="21400602")
                elif len(self.sample_info[key]) < 2:
                    self.set_error('样本%s对应的R1,R2序列不全,请核实！', variables=(key), code="21400603")



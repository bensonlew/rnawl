# -*- coding: utf-8 -*-
# __author__ = 'sj'

from __future__ import division
import os
import re
import shutil
from Bio import SeqIO
from biocluster.agent import Agent
from biocluster.tool import Tool
# from biocluster.core.exceptions import OptionError
from collections import defaultdict

class LengthStatAgent(Agent):
    """
    从fastq或者fastq文件夹里提取样本的信息
    """
    def __init__(self, parent):
        super(LengthStatAgent, self).__init__(parent)
        options = [
            {"name": "length_dir", "type": "infile", "format": "sample.data_dir"},
            {"name": "file_sample_list", "type": "infile", "format": "sequence.info_txt"}
        ]
        self.add_option(options)
        self.step.add_steps("sample_check")
        self.on('start', self.start_sample_check)
        self.on("end", self.end_sample_check)

    def start_sample_check(self):
        self.step.sample_check.start()
        self.step.update()

    def end_sample_check(self):
        self.step.sample_check.finish()
        self.step.update()

    def check_options(self):

        return True

    def set_resource(self):
        self._cpu = 4
        self._memory = "4G"

class LengthStatTool(Tool):
    def __init__(self, config):
        super(LengthStatTool, self).__init__(config)
        self.longest = ""
        self.allowed_step = [20, 50, 100, 200]

    def run(self):
        super(LengthStatTool, self).run()
        self._create_reads_len_info()
        self.end()

    def _create_reads_len_info(self):
        """
        生成4个reads_len_info文件
        """
        tmp_list = os.listdir(self.option("length_dir").prop["path"])
        length_list = [os.path.join(self.option("length_dir").prop["path"], i) for i in tmp_list]
        tmp_dir = os.path.join(self.work_dir, "output", "reads_len_info")
        if not os.path.exists(tmp_dir):
            os.mkdir(tmp_dir)
        # 寻找最长的序列
        self.logger.info("开始寻找所有样本中最长的序列")
        max_list = list()
        with open(self.option("file_sample_list").prop["path"], "r") as r:
            r.readline()
            for line in r:
                line = line.strip()
                lst =line.split("\t")
                max_length = lst[-1]
                max_list.append(max_length)
        self.longest = max(max_list)
        self.logger.info("最长的序列寻找完毕，长度为" + str(self.longest))
        for i in self.allowed_step:
            self._write_head(i)
        for path in length_list:
            self._len_stat(path)

    def _write_head(self, step):
        """
        往reads_len_info文件里输出表头

        :param step:文件的步长
        """
        file_name = os.path.join(self.work_dir, "output", "reads_len_info",
                                 "step_" + str(step) + ".reads_len_info.txt")
        with open(file_name, "w") as f:
            f.write("sample" + "\t")
            col = int(self.longest) + step
            head = list()
            for i in range(step, col, step):
                if step == 1:
                    str_ = str(i - step + 1)
                else:
                    str_ = str(i - step + 1) + "-" + str(i)
                head.append(str_)
            f.write("\t".join(head) + "\n")

    def _len_stat(self, file):
        """
        统计一个fasta文件的长度分布, 往不同的step文件里输出一行(一个样本的长度分布）
        step_1: 用于记录步长为1的reads的分布信息
        step_20: 用于记录步长为20的reads的分布信息 例如区间21-40 对应的应该是step_20[40]
        step_50: 用于记录步长为50的reads的分布信息
        step_100: 用于记录步长为100的reads的分布信息
        """
        sample_name = os.path.basename(file)
        sample_name = re.sub(r"\.length_file$", r"", sample_name)
        self.step_1 = defaultdict(int)
        self.step_20 = defaultdict(int)
        self.step_50 = defaultdict(int)
        self.step_100 = defaultdict(int)
        self.step_200 = defaultdict(int)
        with open(file, "r") as r:
            for line in r:
                line = line.strip()
                len_ = int(line)
                self.logger.info(len_)
                for i in self.allowed_step:
                    self._find_range(len_, i, eval("self.step_" + str(i)))
            for mystep in self.allowed_step:
                self._write_len_info(mystep, eval("self.step_" + str(mystep)), sample_name)

    def _write_len_info(self, step, dict_, sample_name):
        """
        往step_1.reads_len_info;step_20.reads_len_info;step_50.reads_len_info;step_100.reads_len_info
        输出一行

        :param step: 步长
        :param dict_: 步长对应的字典，长度分布数据
        :param sample_name: 样本名称
        """
        file_name = os.path.join(self.work_dir, "output", "reads_len_info",
                                 "step_" + str(step) + ".reads_len_info.txt")
        with open(file_name, "a") as f:
            temp_list = list()
            temp_list.append(sample_name)
            col = int(self.longest) + step
            for i in range(step, col, step):
                temp_list.append(dict_[i])
            f.write("\t".join(str(x) for x in temp_list) + "\n")

    @staticmethod
    def _find_range(len_, step, dict_):
        """
        计算某一个长度序列应该属于哪个区间，并将相应的dict 加1
        例如某条序列 长度len_为32，要计算步长20时，属于哪个区间，则传入参数应当是(32, 20, step_20)
        最后计算可知32 属于21-40的区间，字典step_20[40]应当加1

        :param len_:  序列的长度
        :param step:  步长
        :param dict_: 需要处理的字典
        """
        i = step
        while True:
            if i // len_ >= 1:  # modified by sj on 20161122
                dict_[i] += 1
                break
            i += step



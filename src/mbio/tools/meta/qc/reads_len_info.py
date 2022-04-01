# -*- coding: utf-8 -*-
# __author__ = 'xuting'
from __future__ import division
import math
import os
import re
from collections import defaultdict
from Bio import SeqIO
from biocluster.agent import Agent
from biocluster.tool import Tool
from mbio.files.sequence.fasta import FastaFile
from biocluster.core.exceptions import OptionError


class ReadsLenInfoAgent(Agent):
    """
    version 1.0
    author: xuting
    last_modify: 2015.11.09
    """
    def __init__(self, parent):
        """
        :param longset:所有的fasta中最长的那条序列
        :param allow_step:允许的步长
        """
        super(ReadsLenInfoAgent, self).__init__(parent)
        options = [{"name": "fasta_path", "type": "infile", "format": "sequence.fasta_dir"}]  # 输入的fasta文件夹
        self.add_option(options)
        self.step.add_steps("seq_len_stat")
        self.on('start', self.start_len_stat)
        self.on('end', self.end_len_stat)

    def start_len_stat(self):
        self.step.seq_len_stat.start()
        self.step.update()

    def end_len_stat(self):
        self.step.seq_len_stat.finish()
        self.step.update()

    def check_options(self):
        """
        参数检测
        :return:
        """
        if not self.option("fasta_path").is_set:
            raise OptionError("输入的fasta文件夹路径没有指定")
        return True

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
            ["reads_len_info", "", "样本长度分布信息文件夹"]
        ])
        result_dir.add_regexp_rules([
            ["\.reads_len_info\.txt$", "xls", "样本长度分布信息文件"]
        ])
        super(ReadsLenInfoAgent, self).end()

    def set_resource(self):
        """
        设置所需资源
        """
        self._cpu = 10
        total = 0
        for f in self.option("fasta_path").prop["fasta_fullname"]:
            total += os.path.getsize(f)
        total = total / (1024 * 1024 * 1024)
        total = total * 4
        total = math.ceil(total)
        self._memory = '{}G'.format(int(total))


class ReadsLenInfoTool(Tool):
    """
    """
    def __init__(self, config):
        super(ReadsLenInfoTool, self).__init__(config)
        self._version = 1.0
        self.longest = ""
        self.allowed_step = [20, 50, 100, 200]

    def _create_reads_len_info(self):
        """
        生成4个reads_len_info文件
        """
        tmp_dir = os.path.join(self.work_dir, "output", "reads_len_info")
        if not os.path.exists(tmp_dir):
            os.mkdir(tmp_dir)
        # 寻找最长的序列
        self.logger.info("开始寻找所有样本中最长的序列")
        max_list = list()
        self.option("fasta_path").get_full_info(os.path.join(self.work_dir, "fasta_work_dir"))
        for fasta in self.option("fasta_path").prop["fasta_fullname"]:
            # 获取每一个fasta的全路径
            myfasta = FastaFile()
            myfasta.set_path(fasta)
            myfasta.get_info()
            max_list.append(int(myfasta.prop["longest"]))
        self.longest = max(max_list)
        self.logger.info("最长的序列寻找完毕，长度为" + str(self.longest))
        for i in self.allowed_step:
            self._write_head(i)
        for fasta in self.option("fasta_path").prop["fasta_fullname"]:
            self.logger.info("开始统计 " + str(fasta) + " 的长度分布")
            self._len_stat(fasta)
            self.logger.info(str(fasta) + " 的长度分布统计完毕")

    def _write_head(self, step):
        """
        往reads_len_info文件里输出表头

        :param step:文件的步长
        """
        file_name = os.path.join(self.work_dir, "output", "reads_len_info",
                                 "step_" + str(step) + ".reads_len_info.txt")
        with open(file_name, "w") as f:
            f.write("sample" + "\t")
            col = self.longest + step
            head = list()
            for i in range(step, col, step):
                if step == 1:
                    str_ = str(i - step + 1)
                else:
                    str_ = str(i - step + 1) + "-" + str(i)
                head.append(str_)
            f.write("\t".join(head) + "\n")

    def _len_stat(self, fasta):
        """
        统计一个fasta文件的长度分布, 往不同的step文件里输出一行(一个样本的长度分布）
        step_1: 用于记录步长为1的reads的分布信息
        step_20: 用于记录步长为20的reads的分布信息 例如区间21-40 对应的应该是step_20[40]
        step_50: 用于记录步长为50的reads的分布信息
        step_100: 用于记录步长为100的reads的分布信息
        """
        sample_name = os.path.basename(fasta)
        sample_name = re.sub(r"\.(fasta|fa)$", r"", sample_name)
        self.step_1 = defaultdict(int)
        self.step_20 = defaultdict(int)
        self.step_50 = defaultdict(int)
        self.step_100 = defaultdict(int)
        self.step_200 = defaultdict(int)
        for seq in SeqIO.parse(fasta, "fasta"):
            len_ = len(seq.seq)
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
            col = self.longest + step
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

    def run(self):
        """
        运行
        """
        super(ReadsLenInfoTool, self).run()
        self._create_reads_len_info()
        self.logger.debug("所有样本的长度统计均已完成，程序即将退出")
        self.end()

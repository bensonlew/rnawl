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

class SampleCheckAgent(Agent):
    """
    从fastq或者fastq文件夹里提取样本的信息
    """
    def __init__(self, parent):
        super(SampleCheckAgent, self).__init__(parent)
        options = [
            {"name": "file_sample_list", "type": "infile", "format": "sequence.info_txt"},  # 上游产生的样本信息txt文件
            {"name": "otu_fasta", "type": "outfile", "format": "sequence.fasta"},
            {"name": "raw_sequence", "type" : "infile", "format": "sequence.raw_sequence_txt"}  # 原始序列文件
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

class SampleCheckTool(Tool):
    def __init__(self, config):
        super(SampleCheckTool, self).__init__(config)
        self.longest = ""
        self.allowed_step = [20, 50, 100, 200]

    def cat_fastas_dir(self, sample_list):
        #seq_nu = {}  # 增加序列重编号功能，避免存在重复的序列id导致后续分析出错 add by zhujuan 2018.04.28
        cat_fasta = self.work_dir + "/cat_meta.fasta"
        if os.path.exists(cat_fasta):
            os.remove(cat_fasta)
        os.mknod(cat_fasta)
        # for path in path_list:
        for path in sample_list:
        #     # path += "/output/fa"
        #     for fasta in os.listdir(path):
        #         fasta = os.path.join(path,fasta)
        #         basename = os.path.basename(fasta)
                # self.logger.info(basename)
                # self.logger.info(fasta)
                # break
            basename = os.path.basename(path)
            sample_name = re.search(r'(.+)\.(fasta|fa)$', basename).group(1)
            #if str(sample_name) not in seq_nu.keys():
                #seq_nu[str(sample_name)] = 0
            with open(cat_fasta, "a") as f:
                # for seq in SeqIO.parse(fasta, "fasta"):
                for seq in SeqIO.parse(path, "fasta"):
                    new_id = str(sample_name) + '_' + str(seq.id)
                    #seq_nu[str(sample_name)] += 1
                    #new_id = str(sample_name) + '_' + str(seq_nu[str(sample_name)])
                    f.write('>' + new_id + "\n")
                    f.write(str(seq.seq) + "\n")
        self.option("otu_fasta").set_path(cat_fasta)
        return cat_fasta

    def cat_fastas_dir1(self, sample_list):
        seq_nu = {}
        cat_fasta = self.work_dir + "/cat_meta.fasta"
        if os.path.exists(cat_fasta):
            os.remove(cat_fasta)
        os.mknod(cat_fasta)
        for path in sample_list:
            basename = os.path.basename(path)
            sample_name = re.search(r'(.+)\.(fasta|fa)$', basename).group(1)
            if str(sample_name) not in seq_nu.keys():
                seq_nu[str(sample_name)] = 0
            with open(cat_fasta, "a") as f:
                for seq in SeqIO.parse(path, "fasta"):
                    seq_nu[str(sample_name)] += 1
                    new_id = str(sample_name) + '_' + str(seq_nu[str(sample_name)])
                    f.write('>' + new_id + "\n")
                    f.write(str(seq.seq) + "\n")
        self.option("otu_fasta").set_path(cat_fasta)
        return cat_fasta

    def run(self):
        super(SampleCheckTool, self).run()
        self.option("file_sample_list").get_info()
        # self.cat_fastas_dir(self.option("file_sample_list").workdir_path)
        # self._create_reads_len_info(self.option("file_sample_list").workdir_path)
        self.cat_fastas_dir1(self.option("file_sample_list").sample_path)
        self._create_reads_len_info(self.option("file_sample_list").length_path)
        self.sample_info(self.option("file_sample_list").prop["path"])
        if self.option("raw_sequence").is_set:
            self.sequence_statistics()
        else:
            self.valid_statistics()
        self.end()

    def _create_reads_len_info(self, length_list):
        """
        生成4个reads_len_info文件
        """
        tmp_dir = os.path.join(self.work_dir, "output", "reads_len_info")
        if not os.path.exists(tmp_dir):
            os.mkdir(tmp_dir)
        # 寻找最长的序列
        self.logger.info("开始寻找所有样本中最长的序列")
        max_list = list()
        """
        self.option("fasta_path").get_full_info(os.path.join(self.work_dir, "fasta_work_dir"))
        for fasta in self.option("fasta_path").prop["fasta_fullname"]:
            # 获取每一个fasta的全路径
            myfasta = FastaFile()
            myfasta.set_path(fasta)
            myfasta.get_info()
            max_list.append(int(myfasta.prop["longest"]))
        """
        with open(self.option("file_sample_list").prop["path"],"r") as r:
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
        """
        for fasta in self.option("fasta_path").prop["fasta_fullname"]:
            self.logger.info("开始统计 " + str(fasta) + " 的长度分布")
            self._len_stat(fasta)
            self.logger.info(str(fasta) + " 的长度分布统计完毕")
        """
        # for path in workdir_path:
        #     path = path + "/output/length"
        #     for file in os.listdir(path):
        #         file = os.path.join(path,file)
        #         self._len_stat(file)
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
        """
        for seq in SeqIO.__parse_details(fasta, "fasta"):
            len_ = len(seq.seq)
        """
        with open(file,"r") as r:
            for line in r:
                line = line.strip()
                len_ = int(line)
                # self.logger.info(len_)
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

    def sample_info(self,info_txt):
        dir_path = self.output_dir + "/samples_info"
        if os.path.exists(dir_path):
            shutil.rmtree(dir_path)
        os.mkdir(self.output_dir + "/samples_info")
        with open(dir_path + "/samples_info.txt","w") as w:
            w.write("sample\treads\tbases\tavg\tmin\tmax\n")
            with open(info_txt,"r") as r:
                r.readline()
                for line in r:
                    line = line.strip()
                    lst = line.split("\t")
                    sample = lst[1]
                    if sample in ['NA']:# add by qingchen.zhang@20190924 用于判断样本名称页数名称
                        self.set_error('样本名称不能为NA，请检查样本名称')
                    w.write(lst[1] + "\t" + lst[3] + "\t" +lst[4] + "\t" + lst[5] + "\t" + lst[6] +  "\t" + lst[7] + "\n")

    def sequence_statistics(self):
        self.logger.info("开始统计有效序列信息valid sequence")
        os.link(self.option("raw_sequence").prop["path"],self.output_dir + "/raw_sequence.txt")
        with open(self.option("raw_sequence").prop["path"],"r") as r:
            with open(self.output_dir + "/valid_sequence.txt","w") as w:
                r.readline()
                line = r.readline()
                line = line.strip()
                lst = line.split("\t")
                w.write(lst[0])
        with open(self.output_dir + "/valid_sequence.txt","a") as a:
            with open(self.output_dir + "/samples_info/samples_info.txt", "r") as r:
                r.readline()
                samples = 0
                sequences = 0
                bases = 0
                for line in r:
                    line = line.strip()
                    lst = line.split("\t")
                    sequences += int(lst[1])
                    bases += int(lst[2])
                    samples += 1
                a.write("\t" + str(samples) + "\t" + str(sequences) + "\t" + str(bases) + "\t" + str(bases/sequences))

    def valid_statistics(self):
        self.logger.info("开始统计有效序列信息valid sequence")
        with open(self.output_dir + "/valid_sequence.txt","a") as a:
            with open(self.output_dir + "/samples_info/samples_info.txt", "r") as r:
                r.readline()
                samples = 0
                sequences = 0
                bases = 0
                for line in r:
                    line = line.strip()
                    lst = line.split("\t")
                    sequences += int(lst[1])
                    bases += int(lst[2])
                    samples += 1
                a.write("-" + "\t" + str(samples) + "\t" + str(sequences) + "\t" + str(bases) + "\t" + str(bases/sequences))

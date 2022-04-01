# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'

import os
import re
import shutil
from Bio import SeqIO
from biocluster.agent import Agent
from biocluster.tool import Tool
from collections import defaultdict
from biocluster.core.exceptions import OptionError
from mbio.packages.metaasv.common_function import link_file,link_dir



class SampleMergeAgent(Agent):
    """
    根据上游产生的统计文件合并fq文件，生成新的raw序列信息文件
    """
    def __init__(self, parent):
        super(SampleMergeAgent, self).__init__(parent)
        options = [
            {"name": "info_path", "type": "infile", "format": "sequence.profile_table"},  # 合并后的统计文件
            {"name": "info_temp", "type": "infile", "format": "sequence.profile_table"},  # 未合并的统计文件
            {"name": "raw_sequence", "type" : "infile", "format": "sequence.raw_sequence_txt"}, # 原始序列文件
            # {"name": "in_fasta", "type": "infile", "format": "sequence.fasta_dir"},  # 输入拆分后的样本fasta序列文件
            {"name": "sample_info", "type": "infile", "format": "sequence.profile_table"}, ## 样本检测的结果
            {"name": "otu_fasta", "type": "outfile", "format": "sequence.fasta"},
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("info_path").is_set:
            raise OptionError("必须设置输入的info_path")
        if not self.option("info_temp").is_set:
            raise OptionError("必须设置参数info_temp")
        return True

    def set_resource(self):
        self._cpu = 4
        self._memory = "10G"


class SampleMergeTool(Tool):
    def __init__(self, config):
        super(SampleMergeTool, self).__init__(config)
        self.longest = ""
        self.allowed_step = [20, 50, 100, 200]

    def run(self):
        super(SampleMergeTool, self).run()
        self.get_merge_fasta()
        self.cat_fastas_dir1(self.sample_path_list)
        self._create_reads_len_info(self.length_list)
        self.sample_info(self.option("info_path").prop["path"])
        if self.option("raw_sequence").is_set:
            self.sequence_statistics()
        else:
            self.valid_statistics()
        self.end()

    def get_merge_fasta(self):
        """
        根据样本重命名生成的info_temp.xls和info_path得到合并后的fasta文件
        得到重命名后的length分布
        :return:
        """
        self.logger.info("开始运行合并fasta的样本和合并length")
        self.sample_path_list = []
        self.sample_list = []
        self.length_list = []
        fasta_path = os.path.join(self.work_dir, "fasta_dir")
        if os.path.exists(fasta_path):
            shutil.rmtree(fasta_path)
        os.mkdir(fasta_path)
        length_path = os.path.join(self.work_dir, "length")
        if os.path.exists(length_path):
            shutil.rmtree(length_path)
        os.mkdir(length_path)
        filename_origin_dict = {}
        with open(self.option("sample_info").prop['path'], "r") as fm:
            fm.readline()
            for lin in fm:
                sp_lin = lin.strip().split("\t")
                file_name = os.path.basename(sp_lin[0])
                origin_name = sp_lin[1]
                origin_fa_path = sp_lin[2] + "/output/fa/" + origin_name + ".fasta"
                origin_length_path = sp_lin[2] + "/output/length/" + origin_name + ".length_file"
                if file_name not in filename_origin_dict:
                    filename_origin_dict[file_name] = [{origin_name:(origin_fa_path, origin_length_path)}]
                else:
                    aa_list = filename_origin_dict[file_name]
                    sample_info = {origin_name:(origin_fa_path, origin_length_path)}
                    if sample_info not in aa_list:
                        aa_list.append(sample_info)
                        filename_origin_dict[file_name] = aa_list
        ## 这里为了合并不同样本的样本信息
        with open(self.option("info_temp").prop['path'], 'r') as f:
            f.readline()
            for line in f:
                sp_line = line.strip().split("\t")
                origin_name2 = sp_line[0]
                new_name = sp_line[1]
                new_sample_path = os.path.join(fasta_path, new_name + ".fa")
                new_sample_length_path = os.path.join(length_path, new_name + ".length_file")
                if new_sample_path not in self.sample_path_list:
                    self.sample_path_list.append(new_sample_path)
                if new_sample_length_path not in self.length_list:
                    self.length_list.append(new_sample_length_path)
                if new_name not in self.sample_list:
                    self.sample_list.append(new_name)
                file_name2 = sp_line[2]
                if file_name2 in filename_origin_dict:
                    file_list = filename_origin_dict[file_name2]
                    file_sample= [x.keys()[0] for x in file_list]
                    if origin_name2 in file_sample:
                        origin_file_path = [y.values()[0][0] for y in file_list if y.keys()[0] == origin_name2][0]
                        origin_length_path = [y.values()[0][1] for y in file_list if y.keys()[0] == origin_name2][0]
                        with open(new_sample_path, "a") as f:
                            for seq in SeqIO.parse(origin_file_path, "fasta"):
                                new_id = str(new_name) + '_' + str(seq.id)
                                f.write('>' + new_id + "\n")
                                f.write(str(seq.seq) + "\n")
                        os.system("cat {} >> {}".format(origin_length_path, new_sample_length_path))
        self.logger.info("合并完成fasta的数据！")

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
        with open(self.option("info_path").prop["path"],"r") as r:
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
                    sample = lst[0]
                    if sample in ['NA']:# add by qingchen.zhang@20190924 用于判断样本名称页数名称
                        self.set_error('样本名称不能为NA，请检查样本名称')
                    w.write(lst[0] + "\t" + lst[1] + "\t" +lst[2] + "\t" + lst[3] + "\t" + lst[4] +  "\t" + lst[5] + "\n")

    def sequence_statistics(self):
        self.logger.info("开始统计有效序列信息valid sequence")
        if os.path.exists(self.output_dir + "/raw_sequence.txt"):
            os.remove(self.output_dir + "/raw_sequence.txt")
        os.link(self.option("raw_sequence").prop["path"],self.output_dir + "/raw_sequence.txt")
        with open(self.option("raw_sequence").prop["path"],"r") as r:
            with open(self.output_dir + "/valid_sequence.txt","w") as w:
                r.readline()
                line = r.readline()
                line = line.strip()
                lst = line.split("\t")
                w.write(lst[0])
        with open(self.output_dir + "/valid_sequence.txt","a") as a:
            with open(self.option("info_path").prop["path"], "r") as r:
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
            with open(self.option("info_path").prop["path"], "r") as r:
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

    def end(self):
        super(SampleMergeTool, self).end()
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



class SampleCheckAgent(Agent):
    """
    根据上游产生的统计文件合并fq文件，生成新的raw序列信息文件
    """
    def __init__(self, parent):
        super(SampleCheckAgent, self).__init__(parent)
        options = [
            {"name": "info_path", "type": "infile", "format": "sequence.profile_table"},  # 合并后的统计文件
            {"name": "info_temp", "type": "infile", "format": "sequence.profile_table"},  # 未合并的统计文件
            {"name": "raw_sequence", "type" : "infile", "format": "metaasv.raw_sequence_txt"}, # 原始序列文件
            {"name": "in_fastq", "type": "infile", "format": "sequence.fastq_dir"},  # 输入拆分后的样本fastq序列文件
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("info_path").is_set:
            raise OptionError("必须设置输入的info_path")
        if not self.option("info_temp").is_set:
            raise OptionError("必须设置参数info_temp")
        if not self.option("in_fastq").is_set:
            raise OptionError("必须设置参数in_fastq序列文件夹")
        return True

    def set_resource(self):
        self._cpu = 4
        self._memory = "10G"


class SampleCheckTool(Tool):
    def __init__(self, config):
        super(SampleCheckTool, self).__init__(config)
        self.longest = ""
        self.allowed_step = [20, 50, 100, 200]

    def run(self):
        super(SampleCheckTool, self).run()
        if len(os.listdir(self.output_dir)) == 0:
            self.merge_fastq()
            if self.option("raw_sequence").is_set:
                self.sequence_statistics()
            else:
                self.valid_statistics()
        self.end()

    def merge_fastq(self):
        """
        合并fastq文件
        :return:
        """
        in_fastq = self.option("in_fastq").prop['path']
        info_path = self.option("info_path").prop['path']
        info_temp = self.option("info_temp").prop['path']
        sample_list = []
        all_dirs = os.listdir(in_fastq)
        sample_file_list = {}
        in_fastq_path = "/".join(in_fastq.split("/")[0:-1])
        list_path = os.path.join(in_fastq_path, "list.txt")
        self.logger.info("list_path: {}".format(list_path))
        if os.path.exists(list_path):
            with open(list_path, "r") as f:
                for line in f:
                    sp_line = line.strip().split("\t")
                    new_sample = sp_line[0]
                    origin_name = sp_line[1]
                    file_name = os.path.basename(sp_line[2])
                    sample_split = new_sample.split(".fastq")
                    if len(sample_split) > 1:
                        sample = sample_split[0]
                        origin = origin_name.split(".fastq")[0]
                    else:
                        sample = new_sample.split(".fq")[0]
                        origin = origin_name.split(".fq")[0]
                    if file_name not in sample_file_list:
                        sample_file_list[file_name] = [{origin : sample}]
                    else:
                        origin_dict_list = sample_file_list[file_name]
                        origin_name_list = [x.keys()[0] for x in origin_dict_list]
                        if origin not in origin_name_list:
                            sample_file_list[file_name].append({origin : sample})

        if os.path.exists(os.path.join(self.work_dir, "fq")):
            shutil.rmtree(os.path.join(self.work_dir, "fq"))
        self.logger.info("sample_file_list: {}".format(sample_file_list))
        link_dir(in_fastq, os.path.join(self.work_dir, "fq"))
        new_sample_file = {}
        with open(info_temp, 'r') as f:###下面是要根据修改的样本名称进行合并fq，如果已经合并了文件，这里不再进行合并
            all_temp_sample_list = []
            lines = f.readlines()
            for line in lines[1:]:
                line = line.strip().split("\t")
                new_sample = line[1]
                origin_sample = line[0]
                file_name2 = line[2]
                if new_sample in ['NA']:
                    self.set_error('样本名称不能为NA，请检查样本名称')
                if new_sample not in all_temp_sample_list:
                    all_temp_sample_list.append(new_sample)
                if file_name2 in sample_file_list:
                    split_sample_list = sample_file_list[file_name2]
                    sample_file_name_li = [x.values()[0] for x in split_sample_list if x.keys()[0] == origin_sample]
                    self.logger.info("sample_file_name_li: {}".format(sample_file_name_li))
                    if len(sample_file_name_li) > 0:
                        sample_file_name = sample_file_name_li[0]
                        if new_sample not in new_sample_file:
                            new_sample_file[new_sample] = [sample_file_name]
                        else:
                            new_sample_list = new_sample_file[new_sample]
                            new_sample_list.append(sample_file_name)
                            new_sample_file[new_sample] = new_sample_list
        self.logger.info("new_sample_file: {}".format(new_sample_file))
        fq_dir = os.path.join(self.work_dir, "fq")
        out_fq_dir = os.path.join(self.output_dir, "fq")
        if os.path.exists(out_fq_dir):
            shutil.rmtree(out_fq_dir)
        os.mkdir(out_fq_dir)
        self.logger.info("all_temp_sample_list : {}".format(all_temp_sample_list))
        for new_sample in all_temp_sample_list:
            # file_list = [] ## 存重命名后原始的样本名称
            # for n_specimen in all_new_specimen.keys():
            #     select_specimen = all_new_specimen[n_specimen]
            #     if new_sample in [select_specimen]:
            #         file_list.append(n_specimen)
            file_list = new_sample_file[new_sample]
            if len(file_list) == 1:##没有合并样本的情况
                self.logger.info("file_list : {}".format(file_list[0]))
                origin_path = os.path.join(fq_dir, file_list[0]+ ".fastq")
                new_path = os.path.join(out_fq_dir, new_sample + ".fq")
                origin_path2 = os.path.join(fq_dir, file_list[0]+ ".fq")
                if os.path.exists(origin_path):
                    link_file(origin_path, new_path)
                elif os.path.exists(origin_path2):
                    link_file(origin_path2, new_path)
            elif len(file_list) > 1:##存在合并样本的情况
                new_file_path = os.path.join(out_fq_dir, new_sample + ".fq")
                for file_new in file_list:
                    file_path = os.path.join(fq_dir, file_new + ".fastq")
                    file_path2 = os.path.join(fq_dir, file_new + ".fq")
                    if os.path.exists(file_path):
                        os.system("cat {} >> {}".format(file_path, new_file_path))
                    elif os.path.exists(file_path2):
                        os.system("cat {} >> {}".format(file_path2, new_file_path))

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
        super(SampleCheckTool, self).end()
# -*- coding: utf-8 -*-
# __author__ = 'shijin'

from __future__ import division
import os
import re
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from collections import defaultdict
import resource


class MetaFastqExtractAgent(Agent):
    """
    拆分fq序列
    """
    def __init__(self, parent):
        super(MetaFastqExtractAgent, self).__init__(parent)
        options = [
            {"name": "in_fastq", "type": "infile", "format": "sequence.fastq"},

        ]
        self.add_option(options)
        self.step.add_steps("sample_extract")
        self.on('start', self.start_sample_extract)
        self.on("end", self.end_sample_extract)

    def start_sample_extract(self):
        self.step.sample_extract.start()
        self.step.update()

    def end_sample_extract(self):
        self.step.sample_extract.finish()
        self.step.update()

    def check_options(self):
        if not self.option("in_fastq").is_set:
            raise OptionError("参数in_fastq不能为空")

    def set_resource(self):
        self._cpu = 2
        self._memory = "5G"


class Sample(object):
    def __init__(self, name, fastq_dir='./mytest2/', length_dir='./mytest2/'):
        self.name = name
        self.total_base_num = 0
        self.seqs_len = defaultdict(int)
        self._new_fastq_file = open(fastq_dir + '/' + self.name + '.fastq', 'w')
        self.length_dir = length_dir

    def add_new_fastq(self, seq1, seq2, seq3, seq4):
        self._new_fastq_file.write('{}{}{}{}'.format(seq1, seq2, seq3, seq4))
        length = len(seq2) - 1  # 减去一个回车符
        self.seqs_len[length] += 1
        self.total_base_num += length

    def get_info(self):
        self.seqs_num = sum(self.seqs_len.values())
        self.seqs_mean_len = self.total_base_num / self.seqs_num
        self.seqs_min_len = min(self.seqs_len.keys())
        self.seqs_max_len = max(self.seqs_len.keys())

    def write_len(self):
        self._new_seq_len_file = open(self.length_dir + '/' + self.name + '.length_file', 'w')
        for length in sorted(self.seqs_len.keys()):
            self._new_seq_len_file.write('\n'.join([str(length)] * self.seqs_len[length]) + '\n')

    def close_all(self):
        self._new_fastq_file.close()
        self._new_seq_len_file.close()


class MetaFastqExtractTool(Tool):
    def __init__(self, config):
        super(MetaFastqExtractTool, self).__init__(config)
        self.samples = {}
        self.logger.debug(resource.getrlimit(resource.RLIMIT_NOFILE))
        self.pigz_path = "bioinfo/seq/pigz-2.4/pigz"

    def parse_fastq(self, f_path):
        with open(f_path) as fastq:
            for line in fastq:
                m = re.match("@(.+)_(\d+)", line)
                if not m:
                    self.set_error("fastq文件格式不符合要求，第一行形式应为应为@样本名_序列号")
                    raise Exception('fastq文件格式不符合要求，第一行形式应为应为@样本名_序列号')
                sample_name = m.group(1)
                sample = self.return_sample(sample_name)
                try:
                    sample.add_new_fastq(line, next(fastq), next(fastq), next(fastq))  # 将序列文件拆成多个单样本fq文件
                except:
                    self.set_error("fastq文件缺失，请检查后几行文件是否完整")
                    raise Exception("fastq文件缺失，请检查后几行文件是否完整")
        with open('info.txt', 'w') as info:
            info.write('#file\tsample\tworkdir\tseqs_num\tbase_num\tmean_length\tmin_length\tmax_length\n')
            for sample in self.samples.values():
                sample.get_info()
                sample.write_len()
                sample.close_all()
                info.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(self.option("in_fastq").prop["path"],
                                                                     sample.name,
                                                                     self.work_dir,
                                                                     sample.seqs_num,
                                                                     sample.total_base_num,
                                                                     sample.seqs_mean_len,
                                                                     sample.seqs_min_len,
                                                                     sample.seqs_max_len))

    def return_sample(self, sample_name):
        if sample_name in self.samples:
            return self.samples[sample_name]
        fastq_dir = self.output_dir + "/fastq"
        length_dir = self.output_dir + "/length"
        if not os.path.exists(fastq_dir):
            os.mkdir(fastq_dir)
        if not os.path.exists(length_dir):
            os.mkdir(length_dir)
        sample = Sample(sample_name, fastq_dir=fastq_dir, length_dir=length_dir)
        self.samples[sample_name] = sample
        return sample

    def file_gz(self):
        """
        对文件进行压缩
        """
        file_list = []
        for f in os.listdir(self.output_dir + "/fastq"):
            if f.endswith("gz"):
                os.remove(os.path.join(self.output_dir, "fastq/" + f))
            else:
                file_list.append(os.path.join(self.output_dir, "fastq/" + f))
        cmd = "{} -k {}".format(self.pigz_path, " ".join(file_list))
        command = self.add_command("pigz_file", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("pigz_file压缩文件运行完成")
        else:
            self.set_error("pigz_file压缩文件运行失败")

    def run(self):
        super(MetaFastqExtractTool, self).run()
        self.parse_fastq(self.option("in_fastq").prop["path"])
        self.file_gz()
        self.end()

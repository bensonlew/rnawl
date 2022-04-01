# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'

from __future__ import division
import os
import re
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from collections import defaultdict
import resource


class FastqSampleExtractAgent(Agent):
    """
    从fastq或者fastq文件夹里提取样本的信息，拆分文件
    """
    def __init__(self, parent):
        super(FastqSampleExtractAgent, self).__init__(parent)
        options = [
            {"name": "in_fastq", "type": "infile", "format": "sequence.fastq"},
            {"name": "file_sample_list", "type": "outfile", "format": "sequence.info_txt"},
            {"name": "out_fq", "type": "outfile", "format": "sequence.fastq_dir"},
            {"name": "length_dir", "type": "outfile", "format": "sequence.length_dir"}
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
        self._memory = "10G"


class Sample(object):
    def __init__(self, name, fasta_dir='./mytest2/', length_dir='./mytest2/'):
        self.name = name
        self.total_base_num = 0
        self.seqs_len = defaultdict(int)
        self._new_fasta_file = ""
        self._new_seq_len_file = ""
        self.length_dir = length_dir
        self.fasta_dir = fasta_dir

    def add_new_fasta(self, seq, seq_name):
        if self._new_fasta_file is None or self._new_fasta_file=="" :
            self._new_fasta_file = open(self.fasta_dir + '/' + self.name + '.fasta', 'w')
        self._new_fasta_file.write('>{}\n{}'.format(seq_name, seq))
        length = len(seq) - 1  # 减去一个回车符
        self.seqs_len[length] += 1
        self.total_base_num += length

    def get_info(self):
        self.seqs_num = sum(self.seqs_len.values())
        self.seqs_mean_len = self.total_base_num / self.seqs_num
        self.seqs_min_len = min(self.seqs_len.keys())
        self.seqs_max_len = max(self.seqs_len.keys())

    def write_len(self):
        if self._new_seq_len_file is None or self._new_seq_len_file=="" :
            self._new_seq_len_file = open(self.length_dir + '/' + self.name + '.length_file', 'w')
        for length in sorted(self.seqs_len.keys()):
            self._new_seq_len_file.write('\n'.join([str(length)] * self.seqs_len[length]) + '\n')

    def close_all(self):
        self._new_fasta_file.close()
        self._new_seq_len_file.close()


class FastqSampleExtractTool(Tool):
    def __init__(self, config):
        super(FastqSampleExtractTool, self).__init__(config)
        self.samples = {}
        self.logger.debug(resource.getrlimit(resource.RLIMIT_NOFILE))

    def parse_fastq(self, f_path):
        if self.option("in_fastq").is_gz:  # 增加壓縮文件判斷 by ghd @20180712
            opt_f_path = self.option("in_fastq").unzipfile
        else:
            opt_f_path = self.option("in_fastq").prop['path']
        with open(f_path) as fastq:
            for line in fastq:
                if '.' in line:
                    line = line.replace('.', '_')
                m = re.match("@(.+)_(\d+)", line)
                if not m:
                    self.set_error("fastq文件格式不符合要求，第一行形式应为应为@样本名_序列号,或@样本名.序列号")
                sample_name = m.group(1)
                seq_name = m.group(2)
                sample = self.return_sample(sample_name)
                seq = next(fastq)
                sample.add_new_fasta(seq, seq_name)
                try:
                    _add = next(fastq)
                    _quality = next(fastq)
                    self.split_fastq(sample_name, line,seq, _add, _quality)
                except:
                    self.set_error("fastq文件缺失，请检查后几行文件是否完整")
        with open('info.txt', 'w') as info:
            info.write('#file\tsample\tworkdir\tseqs_num\tbase_num\tmean_length\tmin_length\tmax_length\n')
            for sample in self.samples.values():
                sample.get_info()
                sample.write_len()
                sample.close_all()
                info.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(opt_f_path,
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
        fasta_dir = self.output_dir + "/fa"
        length_dir = self.output_dir + "/length"
        if not os.path.exists(fasta_dir):
            os.mkdir(fasta_dir)
        if not os.path.exists(length_dir):
            os.mkdir(length_dir)
        sample = Sample(sample_name, fasta_dir=fasta_dir, length_dir=length_dir)
        self.samples[sample_name] = sample
        return sample

    def split_fastq(self, sample_name, first, second, third, fourth):
        fastq_dir = self.output_dir + "/fq"
        if not os.path.exists(fastq_dir):
            os.mkdir(fastq_dir)
        fastq_file_path = os.path.join(fastq_dir, sample_name +".fq")
        if not os.path.exists(fastq_file_path):
            self.fastq_file = open(fastq_file_path, 'w')
            self.fastq_file.write("{}{}{}{}".format(first, second, third, fourth))
        else:
            self.fastq_file.write("{}{}{}{}".format(first, second, third, fourth))


    def run(self):
        super(FastqSampleExtractTool, self).run()
        if self.option('in_fastq').is_gz:
            self.option("in_fastq")._prepare(self.work_dir)
            self.parse_fastq(self.option("in_fastq").unzipfile)
        else:
            self.parse_fastq(self.option("in_fastq").prop["path"])
        self.end()

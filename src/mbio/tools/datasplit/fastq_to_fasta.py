# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# last_modified: 20180911

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import re


class FastqToFastaAgent(Agent):
    """
    将fastq文件转换成fasta格式的文件,过滤掉过长过短序列，并进行统计
    工具：fastq_to_fasta.pl、trim_seq.pl
    """
    def __init__(self, parent):
        super(FastqToFastaAgent, self).__init__(parent)
        options = [
            {"name": "fastq", "type": "infile", "format": "sequence.fastq"},
            {"name": "min_length", "type": "string", "default": "18"},  # 最短序列
            {"name": "max_length", "type": "string", "default": "32"},  # 最长序列
        ]
        self.add_option(options)

    def check_options(self):
        """
        检查参数设置
        :return:
        """
        if not self.option('fastq'):
            raise OptionError("必须设置输入的fastq文件")

    def set_resource(self):
        """
        设置所需资源
        :return:
        """
        self._cpu = 10
        self._memory = '5G'


class FastqToFastaTool(Tool):
    def __init__(self, config):
        super(FastqToFastaTool, self).__init__(config)
        self.perl = "program/perl-5.24.0/bin/perl"
        self.fastq_to_fasta_path = self.config.PACKAGE_DIR + "/datasplit/fastq_to_fasta.pl"
        self.trim_seq = self.config.PACKAGE_DIR + "/datasplit/trim_seq.pl"
        self.pigz_path = "bioinfo/seq/pigz-2.4/pigz"

    def run(self):
        super(FastqToFastaTool, self).run()
        self.fastq_to_fasta()
        self.run_trim_seq()
        self.reads_stat()
        self.file_gz()
        self.end()

    def fastq_to_fasta(self):
        """
        运行fastq_to_fasta.pl程序，将单个fastq文件转换成fasta文件
        """
        self.fasta_name = os.path.join(self.work_dir, os.path.basename(self.option("fastq").prop["path"]) + ".fasta")
        cmd = "{} {} -i {} -o {}".format(self.perl, self.fastq_to_fasta_path, self.option("fastq").prop["path"], self.fasta_name)
        fasta_command = self.add_command("fastq_to_fasta", cmd).run()
        self.wait(fasta_command)
        if fasta_command.return_code == 0:
            self.logger.info("fastq_to_fasta运行完成")
        else:
            self.set_error("fastq_to_fasta运行出错!")

    def run_trim_seq(self):
        """
        trim_seq.pl,对转换成的fasta去过长过短序列和统计
        """
        cmd = "{} {} -f {} -o {}".format(self.perl, self.trim_seq, self.fasta_name, os.path.join(self.output_dir, "fasta"))
        cmd += " -min {} -max {}".format(self.option("min_length"), self.option("max_length"))
        command = self.add_command("trim_seq", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("trim_seq运行完成")
        else:
            self.set_error("trim_seq运行出错!")

    def reads_stat(self):
        """
        根据trim_seq.o文件统计<18nt/>32nt/clean的reads
        """
        with open(os.path.join(self.work_dir, "trim_seq.o"), "r") as f, open(os.path.join(self.output_dir, "reads_stat.xls"), "w") as w:
            for line in f:
                if re.match(r"discarded.*less than 18 reads", line):
                    less_reads = line.strip().split()[1]
                if re.match(r"discarded.*bigger than 32 reads", line):
                    big_reads = line.strip().split()[1]
                if re.match(r"remain.*clean reads", line):
                    clean_reads = line.strip().split()[1]
            w.write("<18nt\t>32nt\tclean reads\n")
            w.write(less_reads + "\t" + big_reads + "\t" + clean_reads + "\n")

    def file_gz(self):
        """
        对文件进行压缩
        """
        file_list = []
        for f in os.listdir(self.output_dir):
            if f.endswith("fasta"):
                file_list.append(os.path.join(self.output_dir, f))
        cmd = "{} -k {}".format(self.pigz_path, " ".join(file_list))
        command = self.add_command("pigz_file", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("pigz_file压缩文件运行完成")
        else:
            self.set_error("pigz_file压缩文件运行失败")

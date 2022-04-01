# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# last_modified: 20171130

"""TrimSeq.pl 用于miRNA质控去过长过短序列"""
import os
import glob
from biocluster.config import Config
from biocluster.tool import Tool
from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError
from Bio import SeqIO

class TrimSeqAgent(Agent):
    def __init__(self, parent=None):
        super(TrimSeqAgent, self).__init__(parent)
        options = [
            {"name": "fasta", "type": "infile", "format": "sequence.fasta"},  # 样本fasta文件
            {"name": "min_length", "type": "string", "default": "18"},  # 最短序列
            {"name": "max_length", "type": "string", "default": "32"},  # 最长序列
            {"name": "trim_fasta", "type": "outfile", "format": "small_rna.fasta"},  # 最长序列
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("fasta").is_set:
            raise OptionError("请设置fasta文件")

    def set_resource(self):
        self._cpu = "2"
        self._memory = "10G"


class TrimSeqTool(Tool):
    def __init__(self, config):
        super(TrimSeqTool, self).__init__(config)
        self._version = 1.0
        self.perl = "program/perl-5.24.0/bin/perl"
        self.trim_seq = self.config.SOFTWARE_DIR + "/datasplit/bin/trim_seq.pl"

    def run_trim_seq(self):
        """TrimSeq.pl 用于miRNA质控去过长过短序列"""
        out = os.path.basename(self.option("fasta").prop["path"]).split(".")[0] + "_trimed.fasta"
        cmd = "{} {} -f {} -o {} -min {} -max {}".format(self.perl, self.trim_seq,\
               self.option("fasta").prop["path"], out, self.option("min_length"), self.option("max_length"))
        command = self.add_command("trim_seq", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("运行TrimSeq.pl完成")
            self.option("trim_fasta", os.path.join(self.work_dir, out))
        else:
            self.set_error("运行TrimSeq.pl出错")

    def fasta_stat(self):
        '''
        统计fasta长度分布
        '''
        length_dict = {}
        for seq in SeqIO.parse(self.option("trim_fasta").prop["path"], "fasta"):
            length_dict[seq.id] =int(len(seq.seq))
        #length_dict = self.option("trim_fasta").get_contig_len()
        length_dis = dict()
        for seq_id, seq_len in length_dict.items():
            if seq_len in length_dis:
                length_dis[seq_len] += 1
            else:
                length_dis[seq_len] = 1
        out = os.path.basename(self.option("fasta").prop["path"]).split(".")[0] + "_trimed.length.txt"
        with open(out, 'w') as len_w:
            for seq_len in range(18,33):
                seq_num = 0
                if seq_len in length_dis:
                    seq_num = length_dis[seq_len]
                len_w.write("{}\t{}\n".format(seq_len, seq_num))

    def set_output(self):
        files = glob.glob(r'*_trimed.fasta') + glob.glob(r'*_trimed.length.txt')
        for f in files:
            if os.path.exists(os.path.join(self.output_dir, f)):
                os.remove(os.path.join(self.output_dir, f))
            os.link(os.path.join(self.work_dir, f), os.path.join(self.output_dir, f))
        self.logger.info("Done!")

    def run(self):
        super(TrimSeqTool, self).run()
        self.run_trim_seq()
        self.fasta_stat()
        self.set_output()
        self.end()

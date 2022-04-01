# -*- coding: utf-8 -*-
# __author__ = 'haidong.gu'
# __last_modify__ = '2019/5/6'
from biocluster.agent import Agent
from biocluster.tool import Tool
import os
from biocluster.core.exceptions import OptionError
from Bio import SeqIO


class SplitEndsAgent(Agent):
    """
    DemoAgent:
    version 1.0
    """

    def __init__(self, parent):
        super(SplitEndsAgent, self).__init__(parent)
        options = [
            {"name": "fa", "type": "infile", "format": "sequence.fasta", "required": True},
            {"name": "len", "type": "int", "default": 1000},
            {"name": "filter_len", "type": "int", "default": 1000},  # 少于此长度的数据不进行两端截取, 放到out_short参数中
            {"name": "out_fa", "type": "outfile", "format": "sequence.fasta"},
            {"name": "out_short", "type": "outfile", "format": "sequence.fasta"},
            {"name": "out_leave", "type": "outfile", "format": "sequence.fasta"}
        ]
        self.add_option(options)

    def check_options(self):
        """
        检查参数是否正确
        """
        # modified check_option
        return True

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 2
        self._memory = "10G"


class SplitEndsTool(Tool):
    def __init__(self, config):
        super(SplitEndsTool, self).__init__(config)
        self.leave_seq_list = []  # 太短的序列，直接拿出来
        self.split_seq_list = []  # 满足长度要求的序列
        self.short_seq_list = []  # 不满足长度要求，但高于比对长度的序列

    def run_split(self):
        """
        description
        :return:
        """
        cut_len = self.option("len")
        for seq_record in SeqIO.parse(self.option("fa").prop["path"], "fasta"):
            if len(seq_record) < cut_len:
                self.leave_seq_list.append(seq_record)
            if len(seq_record) < self.option("filter_len") and len(seq_record) >= cut_len:
                self.short_seq_list.append(seq_record)  # 小于最小比对长度阈值的序列直接舍弃掉
            elif len(seq_record) >= self.option("filter_len"):
                seq_record_left = seq_record[:cut_len]
                seq_record_left.id = seq_record.id + "_left"
                seq_record_left.description = ""
                seq_record_right = seq_record[-cut_len:]
                seq_record_right.id = seq_record.id + "_right"
                seq_record_right.description = ""
                self.split_seq_list.append(seq_record_left)
                self.split_seq_list.append(seq_record_right)

    def set_output(self):
        """
        设置输出文件路径
        :return:
        """
        seq_path1 = os.path.join(self.output_dir, "split.fa")
        seq_path2 = os.path.join(self.output_dir, "short.fa")
        seq_path3 = os.path.join(self.output_dir, "leave.fa")
        if len(self.split_seq_list) > 0:
            SeqIO.write(self.split_seq_list, seq_path1, "fasta")
            self.option("out_fa").set_path(seq_path1)
        if len(self.short_seq_list) > 0:
            SeqIO.write(self.short_seq_list, seq_path2, "fasta")
            self.option("out_short").set_path(seq_path2)
        if len(self.leave_seq_list) > 0:
            SeqIO.write(self.leave_seq_list, seq_path3, "fasta")
            self.option("out_leave").set_path(seq_path3)

    def run(self):
        super(SplitEndsTool, self).run()
        self.run_split()
        self.set_output()
        self.end()
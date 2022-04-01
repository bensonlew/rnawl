# -*- coding: utf-8 -*-
# __author__ = 'hao.gao'
# __last_modify__ = '2020/07/13'

from biocluster.agent import Agent
from biocluster.tool import Tool
import os
from biocluster.core.exceptions import OptionError
from Bio import SeqIO


class SplitChrAgent(Agent):
    """
    DemoAgent:主要根据组装序列的是否是线性还是环形，进行分开；如果是环形，不进行下面分析；否则，需要进行下面分析
    version 1.0
    """

    def __init__(self, parent):
        super(SplitChrAgent, self).__init__(parent)
        options = [
            {"name": "fa", "type": "infile", "format": "sequence.fasta", "required": True},
            {"name": "pre_log", "type": "infile", "format": "sequence.profile_table", "required": True},  # 可能没有文件
            {"name": "circle", "type": "outfile", "format": "sequence.fasta"},
            {"name": "line", "type": "outfile", "format": "sequence.fasta"},
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


class SplitChrTool(Tool):
    def __init__(self, config):
        super(SplitChrTool, self).__init__(config)
        self.circle_list = []  # 环形序列
        self.line_list = []  # 线性序列

    def run_split(self):
        """
        description
        :return:
        """
        self.circle = self.get_log(self.option("pre_log").prop['path'])
        for seq_record in SeqIO.parse(self.option("fa").prop["path"], "fasta"):
            if len(seq_record) >= 2000:
                if self.circle[seq_record.id] in ["True", 'true']:
                    self.circle_list.append(seq_record)
                elif self.circle[seq_record.id] in ["False", 'false']:
                    self.line_list.append(seq_record)
            else:
                pass

    def get_log(self, file):
        dict = {}
        with open(file, "r") as f:
            lines = f.readlines()
            for line in lines:
                line = line.strip().split("\t")
                dict[line[1]] = line[3]
        return dict

    def set_output(self):
        """
        设置输出文件路径
        :return:
        """
        seq_path1 = os.path.join(self.output_dir, "circle.fa")
        seq_path2 = os.path.join(self.output_dir, "line.fa")
        if len(self.circle_list) > 0:
            SeqIO.write(self.circle_list, seq_path1, "fasta")
            self.option("circle").set_path(seq_path1)
        if len(self.line_list) > 0:
            SeqIO.write(self.line_list, seq_path2, "fasta")
            self.option("line").set_path(seq_path2)

    def run(self):
        super(SplitChrTool, self).run()
        self.run_split()
        self.set_output()
        self.end()
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


class TrimSeqAgent(Agent):
    def __init__(self, parent=None):
        super(TrimSeqAgent, self).__init__(parent)
        options = [
            {"name": "fasta", "type": "infile", "format": "sequence.fasta"},  # 样本fasta文件
            {"name": "min_length", "type": "string", "default": "18"},  # 最短序列
            {"name": "max_length", "type": "string", "default": "32"},  # 最长序列
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
        self.trim_seq = self.config.SOFTWARE_DIR + "/bioinfo/seq/scripts/trim_seq.pl"

    def run_trim_seq(self):
        """TrimSeq.pl 用于miRNA质控去过长过短序列"""
        out = os.path.basename(self.option("fasta").prop["path"]).split(".")[0] + "_trimed.fasta"
        cmd = "{} {} -f {} -o {} -min {} -max {}".format(self.perl, self.trim_seq,\
               self.option("fasta").prop["path"], out, self.option("min_length"), self.option("max_length"))
        command = self.add_command("trim_seq", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("运行TrimSeq.pl完成")
        else:
            self.set_error("运行TrimSeq.pl出错")

    def set_output(self):
        files = glob.glob(r'*_trimed.fasta')
        for f in files:
            f_ = os.path.join(self.output_dir, f)
            if os.path.exists(f_):
                os.remove(f_)
            os.link(os.path.join(self.work_dir, f), f_)
        self.logger.info("Done!")

    def run(self):
        super(TrimSeqTool, self).run()
        self.run_trim_seq()
        self.set_output()
        self.end()

# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# last_modified: 20171210

"""fastq_extract.pl 从fastq序列中随机抽取n条序列"""
import os
from biocluster.config import Config
from biocluster.tool import Tool
from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError


class SequenceExtractAgent(Agent):
    def __init__(self, parent=None):
        super(SequenceExtractAgent, self).__init__(parent)
        options = [
            {"name": "fastq", "type": "infile", "format": "datasplit.fastq"},  # 样本fastq文件
            {"name": "fastq_r", "type": "infile", "format": "datasplit.fastq"},  # 样本右端fastq文件
            {"name": "fastq_l", "type": "infile", "format": "datasplit.fastq"},  # 样本左端fastq文件
            {"name": "num", "type": "string", "default": "10000"},  # 总共抽取出来的序列数
            {"name": "out_fasta", "type": "outfile", "format": "sequence.fasta"},  # 输出文件fasta
        ]
        self.add_option(options)
        # self.queue = "chaifen"  # 投递到指定的队列chaifen

    def check_options(self):
        if not self.option("fastq").is_set:
            if not self.option("fastq_r").is_set and not self.option("fastq_l").is_set:
                raise OptionError("请设置参数fastq或者fastq_r、fastq_l")

    def set_resource(self):
        self._cpu = "1"
        self._memory = "5G"

    def end(self):
        super(SequenceExtractAgent, self).end()


class SequenceExtractTool(Tool):
    def __init__(self, config):
        super(SequenceExtractTool, self).__init__(config)
        self._version = 1.0
        self.perl = "program/perl-5.24.0/bin/perl"
        self.fastq_extract = self.config.PACKAGE_DIR + "/datasplit/fastq_extract.pl"
        self.python = "program/Python/bin/python"
        self.fastq_extract_fasta = self.config.PACKAGE_DIR + "/datasplit/fastq_extract_fasta.py"

    def run_fastq_extract(self):
        """从fastq中随机抽取出n条序列组成fasta文件"""
        if self.option("fastq").is_set:
            outname = os.path.basename(self.option("fastq").prop["path"]).split(".")[0] + ".fa"
            out = os.path.join(self.work_dir, outname)
            cmd = "{} {} -1 {} -n {} -o {}".format(self.perl, self.fastq_extract, self.option("fastq").prop["path"], self.option("num"), out)
        else:
            outname = os.path.basename(self.option("fastq_l").prop["path"]).split(".")[0] + "_" + os.path.basename(self.option("fastq_r").prop["path"]).split(".")[0] + ".fa"
            out = os.path.join(self.work_dir, outname)
            cmd = "{} {} -1 {} -2 {} -n {} -o {}".format(self.perl, self.fastq_extract, self.option("fastq_l").prop["path"],\
                  self.option("fastq_r").prop["path"], self.option("num"), out)
        command = self.add_command("fastq_extract", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("随机抽取序列fastq_extract运行成功")
        else:
            self.set_error("随机抽取序列fastq_extract运行失败")
        f = os.path.join(self.output_dir, outname)
        if os.path.exists(f):
            os.remove(f)
        os.link(out, f)
        self.option("out_fasta", f)

    def run_fastq_extract_fasta(self):
        if self.option("fastq").is_set:
            outname = os.path.basename(self.option("fastq").prop["path"]).split(".")[0] + ".fa"
            out = os.path.join(self.work_dir, outname)
            cmd = "{} {} -fq1 {} -n {} -out {}".format(self.python, self.fastq_extract_fasta, self.option("fastq").prop["path"],\
                   self.option("num"), out)
        else:
            outname = os.path.basename(self.option("fastq_l").prop["path"]).split(".")[0] + "_" + os.path.basename(self.option("fastq_r").prop["path"]).split(".")[0] + ".fa"
            out = os.path.join(self.work_dir, outname)
            cmd = "{} {} -fq1 {} -fq2 {} -n {} -out {}".format(self.python, self.fastq_extract_fasta, self.option("fastq_l").prop["path"],\
                   self.option("fastq_r").prop["path"], self.option("num"), out)
        command = self.add_command("fastq_extract", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("随机抽取序列fastq_extract运行成功")
        else:
            self.set_error("随机抽取序列fastq_extract运行失败")
        f = os.path.join(self.output_dir, outname)
        if os.path.exists(f):
            os.remove(f)
        os.link(out, f)
        self.option("out_fasta", f)

    def remove_fq_file(self):
        if self.option("fastq").is_set:
            fq_name = os.path.basename(self.option("fastq").prop["path"]).split(".gz")[0]
            fq = os.path.join(self.work_dir, fq_name)
            if os.path.exists(fq):
                os.remove(fq)
        else:
            fq1_name = os.path.basename(self.option("fastq_l").prop["path"]).split(".gz")[0]
            fq2_name = os.path.basename(self.option("fastq_r").prop["path"]).split(".gz")[0]
            fq1 = os.path.join(self.work_dir, fq1_name)
            fq2 = os.path.join(self.work_dir, fq2_name)
            if os.path.exists(fq1):
                os.remove(fq1)
            if os.path.exists(fq2):
                os.remove(fq2)

    def run(self):
        super(SequenceExtractTool, self).run()
        # self.run_fastq_extract()
        self.run_fastq_extract_fasta()
        self.remove_fq_file()
        self.end()

# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# last_modified: 20171210

"""fastq_extract.pl 从fastq序列中随机抽取n条序列"""
import os
from biocluster.config import Config
from biocluster.tool import Tool
from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError


class FastqExtractAgent(Agent):
    def __init__(self, parent=None):
        super(FastqExtractAgent, self).__init__(parent)
        options = [
            {"name": "fastq", "type": "infile", "format": "sequence.fastq"},  # 样本fastq文件
            {"name": "fastq_r", "type": "infile", "format": "sequence.fastq"},  # 样本右端fastq文件
            {"name": "fastq_l", "type": "infile", "format": "sequence.fastq"},  # 样本左端fastq文件
            {"name": "num", "type": "string", "default": "10000"},  # 总共抽取出来的序列数
            {"name": "out_fasta", "type": "outfile", "format": "sequence.fasta"},  # 输出文件fasta
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("fastq").is_set:
            if not self.option("fastq_r").is_set and not self.option("fastq_l").is_set:
                raise OptionError("请设置参数fastq或者fastq_r、fastq_l")

    def set_resource(self):
        self._cpu = "1"
        self._memory = "5G"


class FastqExtractTool(Tool):
    def __init__(self, config):
        super(FastqExtractTool, self).__init__(config)
        self._version = 1.0
        self.perl = "program/perl-5.24.0/bin/perl"
        self.fastq_extract = self.config.SOFTWARE_DIR + "/bioinfo/seq/scripts/fastq_extract.pl"

    def run_fastq_extract(self):
        """从fastq中随机抽取出n条序列组成fasta文件"""
        if self.option("fastq").is_set:
            outname = os.path.basename(self.option("fastq").prop["path"]).split(".")[0] + ".fa"
            out = os.path.join(self.work_dir, outname)
            cmd = "{} {} -1 {} -n {} -o {}".format(self.perl, self.fastq_extract, self.option("fastq").prop["path"], self.option("num"), out)
        else:
            outname = os.path.basename(self.option("fastq_r").prop["path"]).split(".")[0] + "_" + os.path.basename(self.option("fastq_l").prop["path"]).split(".")[0] + ".fa"
            out = os.path.join(self.work_dir, outname)
            cmd = "{} {} -1 {} -2 {} -n {} -o {}".format(self.perl, self.fastq_extract, self.option("fastq_r").prop["path"],\
                  self.option("fastq_l").prop["path"], self.option("num"), out)
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

    def run(self):
        super(FastqExtractTool, self).run()
        self.run_fastq_extract()
        self.end()

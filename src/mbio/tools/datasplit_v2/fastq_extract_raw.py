# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# last_modified: 20171130

"""多样性提取原始raw fastq"""
import os
import re
import json
from biocluster.tool import Tool
from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError


class FastqExtractRawAgent(Agent):
    def __init__(self, parent):
        super(FastqExtractRawAgent, self).__init__(parent)
        options = [
            {"name": "fq1", "type": "infile", "format": "datasplit.fastq"},   # '*.raw.valid.1.fq'
            {"name": "fq2", "type": "infile", "format": "datasplit.fastq"},   # '*.raw.valid.2.fq'
            {"name": "seq2sam", "type": "infile", "format": "datasplit.path"},  # '*.raw.seq2sam.stat'
            {"name": "sample_primer", "type": "infile", "format": "datasplit.path"},  # 文库中样本primer信息
        ]
        self.add_option(options)
        # self.queue = "chaifen"  # 投递到指定的队列chaifen

    def check_options(self):
        if not self.option("fq1").is_set:
            raise OptionError("请设置fq1")
        if not self.option("fq2").is_set:
            raise OptionError("请设置fq2")
        if not self.option("seq2sam").is_set:
            raise OptionError("请设置seq2sam")
        # if not self.option("sample_primer").is_set:
        #     raise OptionError("请设置样本的引物信息表")

    def set_resource(self):
        self._cpu = "1"
        self._memory = "5G"
        size = os.path.getsize(self.option("fq1").prop["path"])
        memory = float(size) / 1024 / 1024 / 1024 * 2
        if memory > 100:
            self._memory = "100G"
        elif memory > 5:
            self._memory = "{}G".format(int(memory))

    def end(self):
        super(FastqExtractRawAgent, self).end()


class FastqExtractRawTool(Tool):
    def __init__(self, config):
        super(FastqExtractRawTool, self).__init__(config)
        self._version = 1.0
        self.python = "program/Python/bin/python"
        self.get_raw_seq = self.config.PACKAGE_DIR + "/datasplit/get_raw_seq.py"

    def run_get_raw(self):
        raw_dir = os.path.join(self.work_dir, "raw")
        if not os.path.exists(raw_dir):
            os.mkdir(raw_dir)
        r1_path = self.option("fq1").prop["path"]
        r2_path = self.option("fq2").prop["path"]
        if r1_path.endswith(".fastq.gz"):
            new_r1_path = os.path.join(self.work_dir, os.path.basename(r1_path).split(".fastq.gz")[0] + ".fastq")
            os.system("gunzip -c {} > {}".format(r1_path, new_r1_path))
            r1_path = new_r1_path
        if r2_path.endswith(".fastq.gz"):
            new_r2_path = os.path.join(self.work_dir, os.path.basename(r2_path).split(".fastq.gz")[0] + ".fastq")
            os.system("gunzip -c {} > {}".format(r2_path, new_r2_path))
            r2_path = new_r2_path
        cmd = "{} {} -r1 {} -r2 {} -s {} -o {}".format(self.python, self.get_raw_seq,
              r1_path, r2_path, self.option("seq2sam").prop["path"], raw_dir)
        command = self.add_command("get_raw_seq", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("get_raw_seq运行成功")
        else:
            self.set_error("get_raw_seq运行失败")

    def file_gz(self):
        """
        对文件进行压缩
        """
        raw_dir = os.path.join(self.work_dir, "raw")
        for f in os.listdir(raw_dir):
            name = f + ".gz"
            if f.endswith(".R1.raw.fastq"):
                sample_name = f.split(".R1.raw.fastq")[0]
                if sample_name in self.primer_info.keys():
                    name = sample_name + "." + self.primer_info[sample_name] + ".R1.raw.fastq.gz"
            elif f.endswith(".R2.raw.fastq"):
                sample_name = f.split(".R2.raw.fastq")[0]
                if sample_name in self.primer_info.keys():
                    name = sample_name + "." + self.primer_info[sample_name] + ".R2.raw.fastq.gz"
            old = os.path.join(raw_dir, f)
            new = os.path.join(self.output_dir, name)
            os.system("gzip -c {} > {}".format(old, new))

    def get_sample_primer(self):
        self.primer_info = {}
        if self.option("sample_primer").is_set:
            f = open(self.option("sample_primer").prop["path"], "rb")
            self.primer_info = json.loads(f.read())

    def run(self):
        super(FastqExtractRawTool, self).run()
        for f in os.listdir(self.output_dir):
            os.remove(os.path.join(self.output_dir, f))
        self.get_sample_primer()
        self.run_get_raw()
        self.file_gz()
        self.end()

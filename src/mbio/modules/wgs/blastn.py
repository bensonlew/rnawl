# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# modified 2018.0516

from biocluster.core.exceptions import OptionError
from biocluster.module import Module
import re
import os


class BlastnModule(Module):
    """
    拆分序列，进行blast，缩短blastn的时间
    """
    def __init__(self, work_id):
        super(BlastnModule, self).__init__(work_id)
        options = [
            {"name": "query_fa", "type": "infile", "format": "sequence.fasta"},  # 要比对的序列
            {"name": "lines", "type": "int", "default": 10000},  # 将fasta序列拆分此行数的多个文件
            {"name": "dbname_nsq", "type": "string"},  # 比对库的路径
            {"name": "outfmt", "type": "string", "default": "6"},  # 输出格式 tab格式是6
            {"name": "num_threads", "type": "int", "default": 8},
            {"name": "evalue", "type": "float", "default": 1e-5},  # 增加参数evalue、num_alignments modified by zengjing 20180508
            {"name": "num_alignments", "type": "int", "default": 5},
            {"name": "sample_id", "type": "string"}  # 用于标记输出文件名字，可以不输入
        ]
        self.add_option(options)
        self.end_times = 0

    def check_options(self):
        if not self.option("query_fa").is_set:
            raise OptionError("请设置要比对的序列query_fa", code="24500301")
        if not self.option("dbname_nsq"):
            raise OptionError("请设置比对库的路径", code="24500302")

    def run_splitfasta(self):
        """
        拆fasta序列
        """
        options = {
            "fasta": self.option("query_fa"),
            "lines": self.option("lines")
        }
        self.split_fasta = self.add_tool("sequence.split_fasta")
        self.split_fasta.set_options(options)
        self.split_fasta.on("end", self.run_blastn)
        self.split_fasta.run()

    def run_blastn(self):
        """
        blastn 进行比对
        """
        self.fasta_files = os.listdir(self.split_fasta.output_dir)
        for f in self.fasta_files:
            query_fa = os.path.join(self.split_fasta.output_dir, f)
            options = {
                "query_fa": query_fa,
                "dbname_nsq": self.option("dbname_nsq"),
                "outfmt": self.option("outfmt"),
                "num_threads": self.option("num_threads"),
                "evalue": self.option("evalue"),
                "num_alignments": self.option("num_alignments"),
                "sample_id": f.split("fasta_")[1]
            }
            self.blastn = self.add_tool("wgs.blast_n")
            self.blastn.set_options(options)
            self.blastn.on("end", self.set_output, f.split("fasta_")[1])
            self.blastn.run()

    def set_output(self, event):
        self.end_times += 1
        obj = event["bind_object"]
        sample_id = self.option("sample_id") if self.option("sample_id") else "blast"
        old = os.path.join(obj.output_dir, event["data"] + ".blast")
        new = os.path.join(self.output_dir, sample_id + ".blast.out.xls")
        if os.path.exists(new):
            new1 = os.path.join(self.output_dir, sample_id + ".blast.out.1.xls")
            os.rename(new, new1)
            os.system("cat {} {} > {}".format(old, new1, new))
            os.remove(new1)
        else:
            os.link(old, new)
        if len(self.fasta_files) == self.end_times:
            self.end()

    def run(self):
        self.run_splitfasta()
        super(BlastnModule, self).run()

    def end(self):
        super(BlastnModule, self).end()

# -*- coding: utf-8 -*-
# __author__ = 'zhaozhigang'
# last_modify:2020.10.27

from biocluster.module import Module
import os
import shutil
from biocluster.core.exceptions import OptionError
from mbio.packages.bacgenome.common import sum_stat
import unittest

class RepeatmaskerModule(Module):
    def __init__(self, work_id):
        super(RepeatmaskerModule, self).__init__(work_id)
        options = [
            {"name": "genome_fa", "type": "infile", "format": "sequence.fasta"},  # 输入序列
            {"name": "sample", "type": "string"}, # 样品名
            {"name": "analysis", "type": "string","default":"umcomplete"},
        ]
        self.add_option(options)
        self.repeatmasker_list = []

    def check_options(self):
        if not self.option("genome_fa").is_set:
            raise OptionError("必须设置参数genome_fa", code="21200301")
        if not self.option("sample"):
            raise OptionError("必须设置样品名称", code="21200302")
        return True

    def run_repeatmasker(self):
        self.genome_path = self.option("genome_fa").prop['path']
        self.repeatmasker = self.add_tool("bacgenome.repeatmasker")
        opts = {
            "input_genome": self.genome_path,
        }
        self.repeatmasker.set_options(opts)
        self.repeatmasker.on("end", self.set_output)
        self.repeatmasker.run()

    def set_output(self):
        sample_dir = os.path.join(self.output_dir, self.option("sample"))
        if os.path.exists(sample_dir):
            shutil.rmtree(sample_dir)
        os.mkdir(sample_dir)
        #for i in os.listdir(self.repeatmasker.output_dir):
        #    if i.split(".")[-1] == "gff":
        #        file_name = i.split(".")[0]
        out_file = self.repeatmasker.output_dir + "/Rep.out"
        if os.path.exists(out_file):
            os.link(out_file, sample_dir + "/" + self.option("sample") + ".out")
        gff_file = self.repeatmasker.output_dir + "/Rep.gff"
        if os.path.exists(gff_file):
            os.link(gff_file, sample_dir + "/" + self.option("sample") + ".gff")
        tbl_file = self.repeatmasker.output_dir + "/Rep.tbl"
        if os.path.exists(tbl_file):
            os.link(tbl_file, sample_dir + "/" + self.option("sample") + ".tbl")
        # 增加判断没有结果文件的情况下，将.fna.out的结果输出到结果目录
        if len(os.listdir(self.repeatmasker.output_dir)) == 0:
            fnaoutfile = self.repeatmasker.work_dir + "/" + self.option("sample") + ".fna.out"
            if os.path.exists(fnaoutfile):
                os.link(fnaoutfile, sample_dir + "/" + self.option("sample") + ".fna.out")
        self.end()

    def run(self):
        super(RepeatmaskerModule, self).run()
        self.run_repeatmasker()

    def end(self):
        super(RepeatmaskerModule, self).end()
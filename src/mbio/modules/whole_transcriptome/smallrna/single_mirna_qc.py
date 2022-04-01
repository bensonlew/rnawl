# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# last_modified: 20171130

import os
import re
from biocluster.config import Config
from biocluster.module import Module
from biocluster.core.exceptions import OptionError
import unittest


class SingleMirnaQcModule(Module):
    """
    miRNA质控的module,用于miRNA质控,单个fastq文件进行质控
    是否去除前三个碱基，去接头，去低值，去未知碱基序列、去过长过短序列
    """
    def __init__(self, work_id):
        super(SingleMirnaQcModule, self).__init__(work_id)
        options = [
            {"name": "fastq", "type": "infile", "format": "sequence.fastq"},  # fastq文件
            {"name": "length", "type": "int", "default": 18},  # 最小序列长度，丢弃比此值更短的序列
            {"name": "adapter", "type": "string", "default": 'AGATCGGAAGAGCACACGTC'},  # 接头序列，如果微量建库，则改为AAAAAA
            {"name": "phred_score", "type": "string", "default": "20"},  # Phred得分（在0到40之间）
            {"name": "length_contain", "type": "int", "default": 51},  # 最小序列长度，丢弃比此值更短的序列
            {"name": "fastq_format", "type": "string", "default": "sanger"},  # Phred得分（在0到40之间）
            {"name": "min_length", "type": "string", "default": "18"},  # 最短序列
            {"name": "max_length", "type": "string", "default": "32"},  # 最长序列
            {"name": "cut_left", "type": "int", "default": 0}, # miRNA是否要切除前3bp，保留前51bp
            {"name": "cut_tail", "type": "int", "default": 0},
            {"name": "out_fastq", "type": "outfile", "format": "sequence.fastq"},
            {"name": "raw_fastq", "type": "outfile", "format": "sequence.fastq"},
            {"name": "sample_name", "type": "string", "default": None},
        ]
        self.add_option(options)
        self.cut_left = self.add_tool("whole_transcriptome.smallrna.cut_left51")

    def finish_update(self, event):
        step = getattr(self.step, event['data'])
        step.finish()
        self.step.update()

    def check_options(self):
        if not self.option("fastq").is_set:
            raise OptionError("必须设置fastq文件")
        if self.option("fastq_format") not in ["sanger", "solexa", "illumina", "Q33", "Q64"]:
            raise OptionError("无法识别该类型的质量体系")

    def run_cut_left(self):
        options = {
            "fastq": self.option("fastq"),
            "cut_left": self.option("cut_left"),
            "length_contain": self.option("length_contain")
        }

        self.cut_left.set_options(options)
        self.step.add_steps("cut_left")
        step = getattr(self.step, "cut_left")
        step.start()
        self.step.update()
        self.cut_left.on('end', self.finish_update, "cut_left")
        self.cut_left.on("end", self.run_fastx_clipper)
        self.cut_left.on("end", self.set_output, "cut_left")
        self.cut_left.run()

    def run_fastx_clipper(self):
        quality = '35'
        if self.option("fastq_format") in ["illumina", "Q64"]:
            quality = '64'
        elif self.option("fastq_format") in ["solexa"]:
            quality = '-5'
        options = {
            "length": self.option("length"),
            "adapter": self.option("adapter"),
            "quality": quality,
            "cut_tail": self.option("cut_tail")
        }
        if self.option("cut_left") > 0:
            options["fastq_s"] = self.cut_left.option("cut_fastq")
        else:
            options["fastq_s"] = self.option("fastq")
        self.fastx_clipper = self.add_tool("whole_transcriptome.smallrna.fastx_clipper")
        self.fastx_clipper.set_options(options)
        self.step.add_steps("fastx_clipper")
        step = getattr(self.step, "fastx_clipper")
        step.start()
        self.step.update()
        self.fastx_clipper.on('end', self.finish_update, "fastx_clipper")
        self.fastx_clipper.on('end', self.set_output, 'fastx_clipper')
        self.fastx_clipper.on("end", self.run_dynamic_trim)
        self.fastx_clipper.run()

    def run_dynamic_trim(self):
        quality = 'sanger'
        if self.option("fastq_format") in ["illumina", "Q64"]:
            quality = 'illumina'
        elif self.option("fastq_format") in ["solexa"]:
            quality = 'solexa'
        options = {
            "fastq": self.fastx_clipper.option("clip_s"),
            "quality": quality,
            "phred_score": self.option("phred_score")
        }
        self.dynamic_trim = self.add_tool("whole_transcriptome.smallrna.dynamic_trim")
        self.dynamic_trim.set_options(options)
        self.step.add_steps("dynamic_trim")
        step = getattr(self.step, "dynamic_trim")
        step.start()
        self.step.update()
        self.dynamic_trim.on('end', self.finish_update, "dynamic_trim")
        self.dynamic_trim.on('end', self.set_output, "dynamic_trim")
        self.dynamic_trim.on("end", self.run_fastq_fasta)
        self.dynamic_trim.run()

    def run_fastq_fasta(self):
        options = {
           "fastq_input": self.dynamic_trim.option("out_fastq")
        }
        self.fastq_fasta = self.add_tool("whole_transcriptome.smallrna.fastq_to_fasta")
        self.fastq_fasta.set_options(options)
        self.step.add_steps("fastq_to_fasta")
        step = getattr(self.step, "fastq_to_fasta")
        step.start()
        self.step.update()
        self.fastq_fasta.on('end', self.finish_update, "fastq_to_fasta")
        self.fastq_fasta.on("end", self.run_trim_seq)
        self.fastq_fasta.run()

    def run_trim_seq(self):
        options = {
            "fasta": os.path.join(self.fastq_fasta.output_dir, "fasta"),
            "min_length": self.option("min_length"),
            "max_length": self.option("max_length")
        }
        self.trim_seq = self.add_tool("whole_transcriptome.smallrna.trim_seq")
        self.trim_seq.set_options(options)
        self.step.add_steps("trim_seq")
        step = getattr(self.step, "trim_seq")
        step.start()
        self.step.update()
        self.trim_seq.on('end', self.finish_update, "trim_seq")
        self.trim_seq.on('end', self.set_output, 'trim_seq')
        self.trim_seq.run()

    def set_output(self, event):
        obj = event["bind_object"]
        if event["data"] == "fastx_clipper":
            pass
            '''
            for f in os.listdir(self.fastx_clipper.output_dir):
                if re.search(r"Sample_QC_stat.xls", f):
                    os.link(self.fastx_clipper.output_dir + "/Sample_QC_stat.xls", os.path.join(self.output_dir, f))
            '''
        if event["data"] == "cut_left":
            fastq_name  =os.path.basename(self.cut_left.option("raw_fastq").prop['path'])
            if os.path.exists(os.path.join(self.output_dir, fastq_name)):
                os.remove(os.path.join(self.output_dir, fastq_name))
            os.link(self.cut_left.option("raw_fastq").prop['path'], os.path.join(self.output_dir, fastq_name))
            self.option("raw_fastq", os.path.join(self.output_dir, fastq_name))

        if event["data"] == "dynamic_trim":
            trim_fq_name = os.path.basename(self.dynamic_trim.option("out_fastq").prop['path'])
            if os.path.exists(os.path.join(self.output_dir, trim_fq_name)):
                os.remove(os.path.join(self.output_dir, trim_fq_name))
            os.link(self.dynamic_trim.option("out_fastq").prop['path'], os.path.join(self.output_dir, trim_fq_name))
            self.option("out_fastq", os.path.join(self.output_dir, trim_fq_name))
        if event["data"] == "trim_seq":
            for f in os.listdir(self.trim_seq.output_dir):
                if os.path.exists(os.path.join(self.output_dir, f)):
                    os.remove(os.path.join(self.output_dir, f))
                os.link(os.path.join(self.trim_seq.output_dir, f), os.path.join(self.output_dir, f))
                # 写统计结果
            if os.path.exists(os.path.join(self.output_dir, "Sample_QC_stat.xls")):
                pass
            else:
                with open(self.fastx_clipper.output_dir + "/Sample_QC_stat.xls", 'r') as stat_in, \
                     open(self.trim_seq.work_dir + "/trim_seq.o", 'r') as stat_trim_in, \
                     open(os.path.join(self.output_dir, "Sample_QC_stat.xls"), 'w') as stat_out:
                    trim_short, trim_long, trim_clean = 0, 0, 0
                    for line in stat_trim_in:
                        if re.match(r"discarded.*less than 18 reads", line):
                            trim_short = line.strip().split()[1]
                        if re.match(r"remain.*clean reads", line):
                            trim_clean = line.strip().split()[1]
                        if re.match(r"discarded.*bigger than 32 reads", line):
                            trim_long = line.strip().split()[1]
                    lines = stat_in.readlines()
                    stat_out.write(lines[0])
                    stats = lines[1].split("\t")
                    if self.option("sample_name"):
                        stats[0] = self.option("sample_name")
                    stats[4] = str(int(stats[4]) + int(trim_short))
                    stats[5] = str(int(stats[5]) + int(trim_long))
                    stats[6] = str(int(stats[6]) + int(trim_clean))
                    stat_out.write("\t".join(stats))
            self.end()

    def run(self):
        super(SingleMirnaQcModule, self).run()
        if self.option("cut_left") > 0 or self.option("length_contain") > 0:
            self.run_cut_left()
        else:
            self.run_fastx_clipper()

    def end(self):
        super(SingleMirnaQcModule, self).end()

class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        import datetime

        test_dir = '/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/whole_transcriptome/demo_human/smallrna/rawdata'
        data = {
            "id": "whole_transcriptome.smallrna_qc" + datetime.datetime.now().strftime('%H-%M-%S'),
            "type": "module",
            "name": "whole_transcriptome.smallrna.single_mirna_qc",
            "instant": False,
            "options": dict(
                fastq = test_dir + "/" + "A1_S12_L001_R1_001.fastq"
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    unittest.main()

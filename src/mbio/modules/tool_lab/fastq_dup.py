# -*- coding: utf-8 -*-
# __author__ = 'fengyitong,fuwenyao,qinjincheng'

import os
import shutil
import unittest

import pandas as pd
from Bio import SeqIO

from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from mbio.files.sequence.file_sample import FileSampleFile

FAKE_FASTQ = """@ST-E00575:252:HJ7HCCCXY:7:1101:8674:1713 1:N:0:CGTACG
GNCCCAACTTGCCATCAAGGATATCTATCTCGGCAACCGCTTCGTTAAATGTCTCTTCGTGGTCAGCTTCAATAGCCAATTCAGATCGGAAGAGCACACGTCTGAACTCCAGTCACCGTACGTAATCTCGTATGCCGTCTTCTGCTTGAAA
+
A#AF-7FFJFJJAJJJJJJJJJJJJJJ-AJFJJJJJJJJJJJJJ<JJJJAJJ7FJJFFFJFFFAFJJJJJAJJJJFFJJJJJJFJFFJJFJJJJFJJJ<FJJJJFJFJJJJJFJJJJJJJJJJJJJJ7FJJJJJJJJ7AJJJJJJJFJJJJ
@ST-E00575:252:HJ7HCCCXY:7:1101:9810:1713 1:N:0:CGTACG
CNGTAATCGTTTGTGGCGTTAGAAATAAAGCCTCAGCCGCCCCGACGACAGAGCCTTCCTTACAAACTTGCCAAAAATAATAAAGATGATTGAAATTGATGTGCGACATTCGCATGTTGTTATCCCCAGATCGGAAGAGCCACACGTCTGA
+
A#AFFJJJJJJJJJJJJJJJJJFJJJJ<JJJJJJJJJJJJAJJJJJJJJJJJJJJJ<-JJJJJJJJJJJJAAJJJJJJJJJJJJJJJJFFJJJJJJJJFFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJFJJJJJJJJJJJA
@ST-E00575:252:HJ7HCCCXY:7:1101:12408:1713 1:N:0:CGTACG
CNCACCAATCATCCTGGACTGGCTCTCAATCTCCATCCTGGAGGTGTCCTTTGTTTCTTCCTGAAACATCCCTTCACTCATCCTAAGCAGTCCCTGAGTCCTTCATCCTGAAGTGGCACCATCCTGATACCGTCCTTTAGATCGGAAGAGC
+
A#AFFFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJFJJFFJFJJJJFJJJJJFFJJJJJ<AJJJJJJJJJJJJJJJJJJJFJFJJJ<AFFJFJJJJJJJJJJJJJJJJJJJFJJJJJJJJJJJJJJJJJJJJJJJJJJJF
"""


class FastqDupModule(Module):
    """
    对hiseq的PE或者SE测序数据做统计，包括GC含量，总reads数目等
    modified by fengyitong at 20180808 -- 增加了用tool运行解压缩的功能
    modified by fuwenyao at 20190506 -- 新增rfam比对；新增dup数据数据展示
    """

    def __init__(self, work_id):
        super(FastqDupModule, self).__init__(work_id)
        options = [
            {"name": "fastq_dir", "type": "infile", "format": "sequence.fastq_dir"},  # fastq文件夹
        ]
        self.add_option(options)
        self.samples = {}
        self.tools = []
        self.gz2fastq = self.add_tool('medical_transcriptome.gzfastq2fastq')
        self.step.add_steps("dup")
        self.if_gz = 0
        self.fake_fastq = FAKE_FASTQ
        self.fq_type = set()

    def check_options(self):
        """
        检查参数
        """
        if not self.option("fastq_dir").is_set:
            raise OptionError("需要传入fastq文件或者文件夹", code="23700601")
        if self.option("fastq_dir").is_set:
            list_path = os.path.join(self.option("fastq_dir").prop["path"], "list.txt")
            if not os.path.exists(list_path):
                raise OptionError("缺少list文件", code="23700602")
            with open(list_path, 'r') as list_r:
                list_info = list_r.read()
                if u'.gz' in list_info:
                    self.if_gz = 1
                    for line in list_info.split('\n'):
                        line = line.strip().split('\t')
                        if u'.gz' in line[0] and not os.path.exists(
                                self.option("fastq_dir").prop["path"] + "/" + ".".join(
                                    line[0].split(".")[:-2]) + ".fastq"):
                            with open(self.option("fastq_dir").prop["path"] + "/" + ".".join(
                                    line[0].split(".")[:-2]) + ".fastq", 'w') as fake:
                                fake.write(self.fake_fastq)
                        if len(line) == 3:
                            self.fq_type.add(line[2])
                else:
                    for line in list_info.split('\n'):
                        line = line.strip().split('\t')
                        if len(line) == 3:
                            self.fq_type.add(line[2])
            self.samples = self.get_list()
            row_num = len(open(list_path, "r").readline().split())
            self.logger.info(row_num)
            if row_num != 3:
                raise OptionError("序列list文件应该包括文件名、样本名和左右端说明三列")
            if len(self.fq_type) == 2:
                for s in sorted(self.samples):
                    if self.samples[s]["l"].split(".")[-1] in ["gz"]:
                        self.samples[s]["l"] = ".".join(self.samples[s]["l"].split(".")[:-2]) + ".fastq"
                    if self.samples[s]["r"].split(".")[-1] in ["gz"]:
                        self.samples[s]["r"] = ".".join(self.samples[s]["r"].split(".")[:-2]) + ".fastq"
            if len(self.fq_type) == 1:
                for s in sorted(self.samples):
                    if isinstance(self.samples[s], dict):
                        if self.samples[s]['s'].split(".")[-1] in ["gz"]:
                            self.samples[s]['s'] = ".".join(self.samples[s]['s'].split(".")[:-2]) + ".fastq"
                    else:
                        if self.samples[s].split(".")[-1] in ["gz"]:
                            self.samples[s] = ".".join(self.samples[s].split(".")[:-2]) + ".fastq"

    def get_list(self):
        list_path = os.path.join(self.option("fastq_dir").prop["path"], "list.txt")
        file_sample = FileSampleFile()
        file_sample.set_path(list_path)
        samples = file_sample.get_list()
        self.logger.info(samples)
        return samples

    def finish_update(self, event):
        step = getattr(self.step, event['data'])
        step.finish()
        self.step.update()

    def dup_finish_update(self):
        self.step.dup.finish()
        self.step.update()

    def run(self):
        super(FastqDupModule, self).run()
        if self.if_gz:
            self.ungzfastq()
        else:
            self.run_tools()

    def ungzfastq(self):
        self.logger.info('需要解压fastq文件，正在解压')
        self.gz2fastq.set_options({'fastq_path': self.option('fastq_dir').prop["path"]})
        self.gz2fastq.on("end", self.run_tools)
        self.gz2fastq.run()

    def run_tools(self):
        self.dup_run()
        self.on_rely(self.tools, self.set_output)
        for dup in self.tools:
            dup.run()
        self.logger.info('{}'.format(self.events))
        for eve in self.events.values():
            self.logger.info('{}'.format(eve.is_start))

    def dup_run(self):
        n = 0
        for sample_name in sorted(self.samples):
            options = {}
            self.logger.info(len(self.fq_type))
            if len(self.fq_type) == 2:
                fq_l = os.path.join(self.option("fastq_dir").prop["path"], self.samples[sample_name]["l"])
                fq_r = os.path.join(self.option("fastq_dir").prop["path"], self.samples[sample_name]["r"])
                options = {'sample_name': sample_name,
                           'fastq_l': fq_l,
                           'fastq_r': fq_r,
                           'fq_type': 'PE'}
            elif len(self.fq_type) == 1:
                fq_s = os.path.join(self.option("fastq_dir").prop["path"], self.samples[sample_name]['s'])
                options = {'sample_name': sample_name,
                           'fastq_s': fq_s,
                           'fq_type': 'SE'}
            dup = self.add_tool('tool_lab.fastq_dup')
            self.step.add_steps('dup_{}'.format(n))
            dup.set_options(options)
            step = getattr(self.step, 'dup_{}'.format(n))
            step.start()
            dup.on("end", self.finish_update, "dup_{}".format(n))
            # dup.run()
            self.tools.append(dup)
            n += 1

    def set_output(self):
        self.logger.info("set output")
        if os.path.isdir(self.output_dir):
            shutil.rmtree(self.output_dir)
        os.mkdir(self.output_dir)
        dup_out = []
        for tool in self.tools:
            out_files = os.listdir(tool.output_dir)
            for f in out_files:
                f_path = os.path.join(tool.output_dir, f)
                if "dup" in f:
                    dup_out.append(f_path)
        with open(self.work_dir + "/dup.xls", "w") as w:
            if len(self.fq_type) == 2:
                w.write("Sample\tDup_R1(%)\tDup_R2(%)\tDup_Pair(%)\n")
            else:
                w.write("Sample\tDup_R1(%)\n")
            for file in dup_out:
                sample_name = os.path.basename(file).split(".")[0]
                sample_name_final = sample_name
                with open(file, "r") as f:
                    f.readline()
                    values = f.next().strip('\n').split('\t')
                    values = [str(round(float(i), 4) * 100) for i in values]
                    w.write("{}\t{}\n".format(sample_name_final, '\t'.join(values)))
        if os.path.exists(self.output_dir + "/dup.xls"):
            os.remove(self.output_dir + "/dup.xls")
        os.link(self.work_dir + "/dup.xls", self.output_dir + "/dup.xls")
        self.logger.info("done")
        self.end()

    def end(self):
        super(FastqDupModule, self).end()

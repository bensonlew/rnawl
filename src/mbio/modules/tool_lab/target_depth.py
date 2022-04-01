# -*- coding: utf-8 -*-
# __author__ = 'zhangyitong'

import os
import shutil
import unittest

import pandas as pd
import glob
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


class TargetDepthModule(Module):
    """
    对hiseq的PE或者SE测序数据做统计，包括GC含量，总reads数目等
    modified by fengyitong at 20180808 -- 增加了用tool运行解压缩的功能
    modified by fuwenyao at 20190506 -- 新增rfam比对；新增dup数据数据展示
    """

    def __init__(self, work_id):
        super(TargetDepthModule, self).__init__(work_id)
        options = [
            {"name": "fastq_dir", "type": "infile", "format": "ref_rna_v2.common_dir"},  # fastq文件夹
            {"name": "seq_file", "type": "infile", "format": "ref_rna_v2.fasta"},
            {"name": "method", 'type': 'string', 'default': 'bwa'},
        ]
        self.add_option(options)
        self.samples = {}
        self.align_tools = []
        self.depth_tools = []
        self.gz2fastq = self.add_tool('medical_transcriptome.gzfastq2fastq')
        self.step.add_steps("align", "depth")
        self.if_gz = 0
        self.fake_fastq = FAKE_FASTQ

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
            row_num = len(open(list_path, "r").readline().split())
            self.logger.info(row_num)
            if row_num != 2:
                raise OptionError("序列list文件应该包括样本名和文件名两列")
            with open(list_path, 'r') as list_r:
                list_info = list_r.read()
                if u'.gz' in list_info:
                    self.if_gz = 1
                    for line in list_info.split('\n')[1:]:
                        line = line.split('\t')
                        if len(line) == 2:
                            fq_files = line[1].split(';')
                            self.get_samples(fq_files, line[0])
                        else:
                            continue
                        for each in fq_files:
                            if u'.gz' in each and not os.path.exists(
                                    self.option("fastq_dir").prop["path"] + "/" + ".".join(
                                        each.split(".")[:-2]) + ".fastq"):
                                with open(self.option("fastq_dir").prop["path"] + "/" + ".".join(
                                        each.split(".")[:-2]) + ".fastq", 'w') as fake:
                                    fake.write(self.fake_fastq)
                else:
                    for line in list_info.split('\n')[1:]:
                        line = line.strip().split('\t')
                        if len(line) == 2:
                            fq_files = line[1].split(';')
                            self.get_samples(fq_files, line[0])
        return True

    def get_samples(self, fq_files, sample):
        if sample not in self.samples.keys():
            self.samples[sample] = dict()
        else:
            self.set_error('Duplicated samples were found.')
        print(self.samples)
        if len(fq_files) == 2:
            if fq_files[0].split(".")[-1] in ['gz']:
                self.samples[sample]["l"] = ".".join(fq_files[0].split(".")[:-2]) + ".fastq"
            else:
                self.samples[sample]['l'] = fq_files[0]
            if fq_files[1].split(".")[-1] in ['gz']:
                self.samples[sample]["r"] = ".".join(fq_files[1].split(".")[:-2]) + ".fastq"
            else:
                self.samples[sample]['r'] = fq_files[1]
            self.samples[sample]['type'] = 'PE'
        else:
            if fq_files[0].split(".")[-1] in ['gz']:
                self.samples[sample]["s"] = ".".join(fq_files[0].split(".")[:-2]) + ".fastq"
            else:
                self.samples[sample]['s'] = fq_files[0]
            self.samples[sample]['type'] = 'SE'

    def finish_update(self, event):
        step = getattr(self.step, event['data'])
        step.finish()
        self.step.update()

    def run(self):
        super(TargetDepthModule, self).run()
        if self.if_gz:
            self.ungzfastq()
        else:
            self.run_align()

    def ungzfastq(self):
        self.logger.info('需要解压fastq文件，正在解压')
        self.gz2fastq.set_options({'fastq_path': self.option('fastq_dir').prop["path"]})
        self.gz2fastq.on("end", self.run_align)
        self.gz2fastq.run()

    def run_align(self):
        self.seq_align()
        self.on_rely(self.align_tools, self.run_depth)
        for align in self.align_tools:
            align.run()
        self.logger.info('{}'.format(self.events))
        for eve in self.events.values():
            self.logger.info('{}'.format(eve.is_start))

    def run_depth(self):
        self.samtools_run()
        self.on_rely(self.depth_tools, self.set_output)
        for depth in self.depth_tools:
            depth.run()
        self.logger.info('{}'.format(self.events))
        for eve in self.events.values():
            self.logger.info('{}'.format(eve.is_start))

    def seq_align(self):
        n = 0
        for sample_name in sorted(self.samples):
            options = {}
            if self.samples[sample_name]['type'] == 'PE':
                fq_l = os.path.join(self.option("fastq_dir").prop["path"], self.samples[sample_name]["l"])
                fq_r = os.path.join(self.option("fastq_dir").prop["path"], self.samples[sample_name]["r"])
                options = {'sample_name': sample_name,
                           'fastq_l': fq_l,
                           'fastq_r': fq_r,
                           'seq_file': self.option('seq_file'),
                           'out': self.output_dir}
            elif self.samples[sample_name]['type'] == 'SE':
                fq_s = os.path.join(self.option("fastq_dir").prop["path"], self.samples[sample_name]['s'])
                options = {'sample_name': sample_name,
                           'fastq_s': fq_s,
                           'seq_file': self.option('seq_file'),
                           'out': self.output_dir}
            align = self.add_tool('tool_lab.target_depth.{}'.format(self.option('method')))
            self.step.add_steps('align_{}'.format(n))
            align.set_options(options)
            step = getattr(self.step, 'align_{}'.format(n))
            step.start()
            align.on("end", self.finish_update, "align_{}".format(n))
            self.align_tools.append(align)
            n += 1

    def samtools_run(self):
        n = 0
        for sample_name in sorted(self.samples):
            sam_f = os.path.join(self.output_dir, '{}_map.sam'.format(sample_name))
            options = {'sample_name': sample_name,
                       'in_sam': sam_f,
                       'fq_type': self.samples[sample_name]['type'],
                       'out': self.output_dir}
            depth = self.add_tool('tool_lab.target_depth.samtools_depth')
            self.step.add_steps('depth_{}'.format(n))
            depth.set_options(options)
            step = getattr(self.step, 'depth_{}'.format(n))
            step.start()
            depth.on("end", self.finish_update, "depth_{}".format(n))
            self.depth_tools.append(depth)
            n += 1

    def set_output(self):
        self.logger.info('start set_output at {}'.format(self.__class__.__name__))
        outpath = os.path.join(self.output_dir, 'uploads')
        if os.path.exists(outpath):
            shutil.rmtree(outpath)
        os.mkdir(outpath)
        for tool in self.depth_tools:
            for source in glob.glob(os.path.join(tool.output_dir, '*')):
                link_name = os.path.join(outpath, os.path.basename(source))
                if os.path.exists(link_name):
                    os.remove(link_name)
                os.link(source, link_name)
                self.logger.info('succeed in linking {} to {}'.format(source, link_name))
        self.logger.info('finish set_output at {}'.format(self.__class__.__name__))
        self.end()

    def end(self):
        super(TargetDepthModule, self).end()

#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import shutil
import glob
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from mbio.files.sequence.file_sample import FileSampleFile


class BwaSamtoolsModule(Module):
    """
    denovoRNA比对后质量评估
    version 1.0
    author: qindanhua
    last_modify: 2016.07.13
    """
    def __init__(self, work_id):
        super(BwaSamtoolsModule, self).__init__(work_id)
        options = [
            {"name": "ref_fasta", "type": "infile", "format": "sequence.fasta"},  # denovo时是基因的参考序列
            {"name": "fq_type", "type": "string", "default": ""},  # fq类型，必传
            {"name": "fastq_r", "type": "infile", "format": "sequence.fastq"},  # 右端序列文件
            {"name": "fastq_l", "type": "infile", "format": "sequence.fastq"},  # 左端序列文件
            {"name": "fastq_s", "type": "infile", "format": "sequence.fastq"},  # SE序列文件
            {"name": "fastq_dir", "type": "infile", "format": "sequence.fastq_dir"},  # fastq文件夹
            # {"name": "head", "type": "string", "default": None},  # 设置结果头文件
            # {"name": "sam", "type": "outfile", "format": "align.bwa.sam_dir"},
            # {"name": "sam", "type": "infile", "format": "align.bwa.sam"},    # sam格式文件
            {"name": "out_bam", "type": "outfile", "format": "align.bwa.bam_dir"},  # bam格式输入文件
            {"name": "method", "type": "string", "default": ""}     # samtool工具
        ]
        self.add_option(options)
        self.samples = {}
        self.samtools = []
        self.bwa_tools = []
        self.bwa = None
        self.end_times = 1
        self.step.add_steps('bwa')
        self.ref_link = ""
        self.index = self.add_tool('align.bwa.bwa')

    def check_options(self):
        """
        检查参数
        """
        if not self.option("ref_fasta").is_set:
            raise OptionError("请传入参考序列")
        if self.option("fastq_dir").is_set:
            self.samples = self.get_list()
            list_path = os.path.join(self.option("fastq_dir").prop["path"], "list.txt")
            row_num = len(open(list_path, "r").readline().split())
            self.logger.info(row_num)
            if self.option('fq_type') == "PE" and row_num != 3:
                raise OptionError("PE序列list文件应该包括文件名、样本名和左右端说明三列")
            elif self.option('fq_type') == "SE" and row_num != 2:
                raise OptionError("SE序列list文件应该包括文件名、样本名两列")
        if self.option('fq_type') not in ['PE', 'SE']:
            raise OptionError("请说明序列类型，PE or SE?")
        if not self.option("fastq_dir").is_set and self.option('fq_type') in ["PE"]:
            if not self.option("fastq_r").is_set:
                raise OptionError("请传入PE右端序列文件")
            if not self.option("fastq_l").is_set:
                raise OptionError("请传入PE左端序列文件")
        if not self.option("fastq_dir").is_set and self.option('fq_type') == "SE":
            if not self.option("fastq_s").is_set:
                raise OptionError("请传入SE序列文件")
        return True

    def bwa_finish_update(self):
        self.step.bwa.finish()
        self.step.update()

    def index_finish_update(self):
        self.step.index.finish()
        self.step.update()

    def finish_update(self, event):
        step = getattr(self.step, event['data'])
        step.finish()
        self.step.update()

    def samtools_finish_update(self, event):
        step = getattr(self.step, event['data'])
        step.finish()
        self.step.update()

    def link_ref(self):
        ref_fasta = self.option('ref_fasta').prop["path"]
        self.ref_link = self.work_dir + "/" + os.path.basename(ref_fasta)
        self.logger.info(self.ref_link)
        if os.path.exists(self.ref_link):
            os.remove(self.ref_link)
        os.link(ref_fasta, self.ref_link)

    def index_run(self):
        self.link_ref()
        self.step.add_steps("index")
        self.index.set_options({
            "ref_fasta": self.ref_link,
            "method": "index"
        })
        self.step.index.start()
        self.index.on("end", self.index_finish_update)
        self.index.on("end", self.multi_bwa_run)
        self.index.run()

    def multi_bwa_run(self):
        # self.samples = self.get_list()
        n = 0
        if self.option("fq_type") == "PE":
            for f in self.samples:
                fq_l = os.path.join(self.option('fastq_dir').prop["path"], self.samples[f]["l"])
                fq_r = os.path.join(self.option('fastq_dir').prop["path"], self.samples[f]["r"])
                bwa = self.add_tool('align.bwa.bwa')
                self.step.add_steps('bwa_{}'.format(n))
                bwa.set_options({
                    "ref_fasta": self.ref_link,
                    'fastq_l': fq_l,
                    'fastq_r': fq_r,
                    'fq_type': self.option('fq_type')
                })
                step = getattr(self.step, 'bwa_{}'.format(n))
                step.start()
                bwa.on("end", self.finish_update, "bwa_{}".format(n))
                # bwa.on("end", self.rename, f)
                bwa.on("end", self.samtools_run, f)
                # bwa.run()
                self.bwa_tools.append(bwa)
        elif self.option("fq_type") == "SE":
            for f in self.samples:
                fq_s = os.path.join(self.option('fastq_dir').prop["path"], self.samples[f])
                bwa = self.add_tool('align.bwa.bwa')
                self.step.add_steps('bwa_{}'.format(n))
                bwa.set_options({
                    "ref_fasta": self.ref_link,
                    'fastq_s': fq_s,
                    'fq_type': self.option('fq_type')
                })
                step = getattr(self.step, 'bwa_{}'.format(n))
                step.start()
                bwa.on("end", self.finish_update, "bwa_{}".format(n))
                # bwa.on("end", self.rename, f)
                bwa.on("end", self.samtools_run, f)
                # bwa.run()
                self.bwa_tools.append(bwa)
        if len(self.bwa_tools) == 1:
            # self.bwa_tools[0].on("end", self.multi_samtools_run)
            self.bwa_tools[0].run()
        else:
            for tool in self.bwa_tools:
                tool.run()
            # self.on_rely(self.bwa_tools, self.multi_samtools_run)

    def bwa_single_run(self):
        self.link_ref()
        self.bwa = self.add_tool('align.bwa.bwa')
        if self.option("fq_type") == "PE":
            self.bwa.set_options({
                "ref_fasta": self.ref_link,
                'fastq_l': self.option('fastq_l').prop["path"],
                'fastq_r': self.option('fastq_r').prop["path"],
                'fq_type': self.option('fq_type')
                })
        elif self.option("fq_type") == "SE":
            self.bwa.set_options({
                "ref_fasta": self.ref_link,
                'fastq_s': self.option('fastq_s').prop["path"],
                'fq_type': self.option('fq_type')
                })
        self.step.bwa.start()
        self.bwa.on("end", self.bwa_finish_update)
        self.bwa.on("end", self.samtools_run, 1)
        self.bwa.run()

    def samtools_run(self, event):
        obj = event["bind_object"]
        bwa_output = os.listdir(obj.output_dir)
        old_name = os.path.join(obj.output_dir, bwa_output[0])
        self.logger.info(event["data"] + ".sam")
        f_path = os.path.join(obj.output_dir, event["data"] + ".sam")
        os.rename(old_name, f_path)
        self.logger.info(f_path)
        self.logger.info(self.ref_link)
        samtools = self.add_tool('denovo_rna.gene_structure.samtools')
        self.step.add_steps('samtools_{}'.format(event["data"]))
        samtools.set_options({
            "ref_fasta": self.ref_link,
            "sam": f_path,
            "method": "sort"
        })
        step = getattr(self.step, 'samtools_{}'.format(event["data"]))
        step.start()
        samtools.on("end", self.finish_update, 'samtools_{}'.format(event["data"]))
        samtools.on("end", self.set_output)
        self.logger.info("samRunnnnn")
        samtools.run()
        self.samtools.append(samtools)

    def rename(self, event):
        obj = event["bind_object"]
        for f in os.listdir(obj.output_dir):
            old_name = os.path.join(obj.output_dir, f)
            new_name = os.path.join(obj.output_dir, event["data"] + ".sam")
            os.rename(old_name, new_name)

    def get_list(self):
        list_path = os.path.join(self.option("fastq_dir").prop["path"], "list.txt")
        file_sample = FileSampleFile()
        file_sample.set_path(list_path)
        samples = file_sample.get_list()
        return samples

    def set_output(self):
        self.logger.info("set output")
        if self.end_times < len(self.samples):
            self.end_times += 1
        elif self.end_times == len(self.samples):
            for f in os.listdir(self.output_dir):
                f_path = os.path.join(self.output_dir, f)
                if os.path.isdir(f_path):
                    shutil.rmtree(f_path)
                else:
                    os.remove(f_path)
            bam_dir = os.path.join(self.output_dir, "sorted_bam")
            os.makedirs(bam_dir)
            for sam in self.samtools:
                outfiles = os.listdir(sam.output_dir)
                # self.logger.info(outfiles)
                for f in outfiles:
                    f_path = os.path.join(sam.output_dir, f)
                    target = bam_dir + "/" + f
                    if "sorted.bam" in f:
                        os.link(f_path, target)
            self.option('out_bam').set_path(bam_dir)
            self.logger.info("set output done")
            self.end()

    def run(self):
        if self.option("fastq_dir").is_set:
            self.index_run()
        else:
            self.bwa_single_run()
        super(BwaSamtoolsModule, self).run()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [r".", "", "结果输出目录"],
            [r"./sorted_bam/", "文件夹", "排序后的bam格式比对结果文件输出目录"]
        ])
        super(BwaSamtoolsModule, self).end()

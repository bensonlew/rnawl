#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import re
import shutil
import time
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from collections import defaultdict


class MetagenomicQcModule(Module):
    """
    宏基因组质控模块，主要调用seqprep、sickle软件做质量剪切与去接头
    version 1.0
    author: wangzhoayue
    last_modify: 2017.12.14
    """

    def __init__(self, work_id):
        super(MetagenomicQcModule, self).__init__(work_id)
        options = [
            {"name": "sample_path", "type": "infile", "format": "sequence.file_sample"},
            # fastq路径list.txt文件，第一列路径，第二列样本名，第三列序列类型 l or r
            # {"name": "fq_type", "type": "string", "default": "PE"},  # PE OR SE
            {"name": "clip_dir", "type": "outfile", "format": "sequence.fastq_dir"},  # SE去接头输出结果文件夹
            {"name": "sickle_dir", "type": "outfile", "format": "sequence.fastq_dir"},  # 质量剪切输出结果文件夹(包括左右段)
            {"name": "sickle_r_dir", "type": "outfile", "format": "sequence.fastq_dir"},  # 质量剪切右端输出结果文件夹
            {"name": "sickle_l_dir", "type": "outfile", "format": "sequence.fastq_dir"},  # 质量剪切左端输出结果文件夹
            {"name": "seqprep_dir", "type": "outfile", "format": "sequence.fastq_dir"},  # PE的去接头输出结果文件
            {"name": "fq_s", "type": "outfile", "format": "sequence.fastq"},  # SE所有样本cat集合
            {"name": "fq_r", "type": "outfile", "format": "sequence.fastq"},  # PE所有右端序列样本cat集合
            {"name": "fq_l", "type": "outfile", "format": "sequence.fastq"},  # PE所有左端序列样本cat集合
            {"name": "quality_q", "type": "int", "default": 20},  # 质量剪切碱基质量 sickle
            {"name": "length_q", "type": "int", "default": 50},  # 质量剪切碱基长度  sickle
            {"name": "seq_quality", "type": "int", "default": 20},  # SeqPrep的参数-q
            {"name": "seq_length", "type": "int", "default": 50},  # SeqPrep的参数-L
            {"name": "min_length", "type": "string", "default": "50"},  # 删除短于此值的reads
        ]
        self.add_option(options)
        self.sample = defaultdict(list)  # 存放样本
        self.seqprep = []
        self.clipper = []
        self.sickle = []
        self.end_times = 0
        self.adapt_rate = []
        self.single_sickle_out = []

    def check_options(self):
        """
        检查参数
        """
        if not self.option("sample_path").is_set:
            raise OptionError("需要传入fastq文件或者文件夹")
        row_num = len(open(self.option("sample_path").prop['path'], "r").readline().split())
        if row_num != 3:
            raise OptionError("list文件应该是两列或三列，PE序列应该包括文件名、样本名和左右端说明三列")

    def get_list(self):
        """
        将输入的信息进行分类,{样本：[R1path,R2path]}
        :return:
        """
        with open(self.option('sample_path').prop['path'])as fr:
            for line in fr:
                tmp = line.strip().split('\t')
                if tmp[1] in self.sample.keys():
                    if tmp[2] == 'l':
                        self.sample[tmp[1]].insert(0, tmp[0])
                    else:
                        self.sample[tmp[1]].append(tmp[0])
                else:
                    self.sample[tmp[1]].append(tmp[0])
            for key in self.sample.keys():
                if len(self.sample[key]) > 2:
                    raise Exception('需要质控的序列PE样本{}有重名，请改样本名或分开质控！'.format(key))

    def run(self):
        super(MetagenomicQcModule, self).run()
        self.get_list()
        time.sleep(2)
        self.seqprep_run()

    def seqprep_run(self):
        n = 1
        for f in self.sample:
            seqprep = self.add_tool('datasplit.seq_prep')
            seqprep.set_options({
                "fastq_l": self.sample[f][0],
                "fastq_r": self.sample[f][1],
                "quality": self.option("seq_quality"),
                "length": self.option("seq_length")
            })
            seqprep.on("end", self.adapt, f)
            seqprep.on("end", self.sickle_pe_run, f)
            n += 1
            self.seqprep.append(seqprep)
        if len(self.seqprep) == 1:
            self.seqprep[0].on("end", self.adapt_write)
            self.seqprep[0].run()
        else:
            self.on_rely(self.seqprep, self.adapt_write)
            for tool in self.seqprep:
                tool.run()

    def sickle_pe_run(self, event):
        obj = event["bind_object"]
        self.tool_sickle = self.add_tool('datasplit.sickle_mg')
        self.tool_sickle.set_options({
            "fq_type": 'PE',
            "fastq_l": obj.option("seqprep_l"),  # modified by shijin on 20170623，减少阻塞
            "fastq_r": obj.option("seqprep_r"),
            "quality": self.option("quality_q"),
            "length": self.option("length_q")
        })
        self.tool_sickle.on("end", self.new_set_output, "sickle_" + event["data"])
        self.tool_sickle.on("end", self.remove_sort_reads, event["data"])
        self.tool_sickle.run()
        self.sickle.append(self.tool_sickle)

    def remove_sort_reads(self, event):
        obj = event["bind_object"]
        remove_sort_reads = self.add_tool('datasplit.remove_sort_reads')
        remove_sort_reads.set_options({
            "fastq_r": obj.option("sickle_l"),
            "fastq_l": obj.option("sickle_l"),
            "min_length": self.option("min_length"),
            "sample_name": event["data"]
        })
        remove_sort_reads.on("end", self.new_set_output, "clean_" + event["data"])
        remove_sort_reads.run()

    def adapt(self, event):
        obj = event["bind_object"]
        adapt_file = obj.work_dir + "/adapter.xls"
        if os.path.exists(adapt_file):
            with open(adapt_file, "r") as f:
                f.readline()
                adapt_rate = f.next().split()[-1]
                self.adapt_rate.append(["{}".format(event["data"]), adapt_rate])

    def adapt_write(self):
        with open(self.output_dir + "/adapter.xls", "w") as w:
            for a in self.adapt_rate:
                w.write("{}\t{}\n".format(a[0], a[1]))

    def set_output(self, event):
        obj = event["bind_object"]
        if self.end_times < len(self.sample):
            self.end_times += 1
        for f in os.listdir(obj.output_dir):
            old_name = os.path.join(obj.output_dir, f)
            new_name = os.path.join(obj.output_dir, event["data"] + "_" + f)
            os.rename(old_name, new_name)
        if self.end_times == len(self.sample):
            sickle_dir = os.path.join(self.output_dir, "sickle_dir")
            sickle_r_dir = os.path.join(self.work_dir, "sickle_r_forRSEM")
            sickle_l_dir = os.path.join(self.work_dir, "sickle_l_forRSEM")
            seqprep_dir = os.path.join(self.work_dir, "seqprep_dir")
            clip_dir = os.path.join(self.work_dir, "clip_dir")
            dir_list = [sickle_dir, seqprep_dir, clip_dir, sickle_r_dir, sickle_l_dir]
            for d in dir_list:
                if os.path.exists(d):
                    shutil.rmtree(d)
                os.mkdir(d)
            sickle_out = []
            seqprep_out = []
            clip_out = []
            for sic in self.sickle:
                for f in os.listdir(sic.output_dir):
                    f_path = os.path.join(sic.output_dir, f)
                    sickle_out.append(f_path)
            for seq in self.seqprep:
                for f in os.listdir(seq.output_dir):
                    f_path = os.path.join(seq.output_dir, f)
                    seqprep_out.append(f_path)
            for clip in self.clipper:
                for f in os.listdir(clip.output_dir):
                    f_path = os.path.join(clip.output_dir, f)
                    clip_out.append(f_path)
            with open(os.path.join(sickle_dir, "list.txt"), "w") as w:
                for f in sickle_out:
                    f_name = f.split("/")[-1]
                    if "sickle_r.fastq" in f:
                        sample_name = f_name.split("_sickle_r.fastq")[0]
                        w.write("{}\t{}\t{}\n".format(f_name, sample_name, "r"))
                        os.link(f, os.path.join(sickle_r_dir, f_name))
                    elif "sickle_l.fastq" in f:
                        sample_name = f_name.split("_sickle_l.fastq")[0]
                        w.write("{}\t{}\t{}\n".format(f_name, sample_name, "l"))
                        os.link(f, os.path.join(sickle_l_dir, f_name))
                    elif "sickle_s.fastq" in f:
                        sample_name = f_name.split("_sickle_s.fastq")[0]
                        w.write("{}\t{}\t{}\n".format(f_name, sample_name, "s"))
                    else:
                        self.logger.info(sickle_dir)
                        self.logger.info(os.path.join(sickle_out, f))
                    target_path = os.path.join(sickle_dir, f_name)
                    if os.path.exists(target_path):
                        os.remove(target_path)
                    os.link(f, target_path)
            self.option("sickle_dir", sickle_dir)
            self.end()

    def link_file(self, old, new):
        if os.path.exists(new):
            os.remove(new)
        os.link(old, new)

    def new_set_output(self, event):
        obj = event["bind_object"]
        list_file = os.path.join(self.output_dir, "list.txt")
        if event["data"].startswith("clean_"):
            self.end_times += 1
            with open(list_file, "w+") as w:
                for f in os.listdir(obj.output_dir):
                    old = os.path.join(obj.output_dir, f)
                    if f.endswith(".1.fq.gz"):
                        new = os.path.join(self.output_dir, event["data"].split("clean_")[1] + ".clean.1.fastq.gz")
                        self.link_file(old, new)
                        w.write(event["data"].split("clean_")[1] + ".clean.1.fastq.gz" + "\t" + event["data"].split("clean_")[1] + "\t" + "l" + "\n")
                    elif f.endswith(".2.fq.gz"):
                        new = os.path.join(self.output_dir, event["data"].split("clean_")[1] + ".clean.2.fastq.gz")
                        self.link_file(old, new)
                        w.write(event["data"].split("clean_")[1] + ".clean.2.fastq.gz" + "\t" + event["data"].split("clean_")[1] + "\t" + "r" + "\n")
        else:
            with open(list_file, "w+") as w:
                for f in os.listdir(obj.output_dir):
                    if re.search("sickle_s.fastq.gz", f):
                        old = os.path.join(obj.output_dir, f)
                        new = os.path.join(self.output_dir, event["data"].split("sickle_")[1] + ".s.clena.fq.gz")
                        self.link_file(old, new)
                        w.write(event["data"].split("sickle_")[1] + ".s.clena.fq.gz" + "\t" + event["data"].split("sickle_")[1] + "\t" + "s" + "\n")
        if self.end_times == len(self.sample):
            self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        # result_dir.add_relpath_rules([
        #     [r".", "", "结果输出目录"],
        #     [r"./seqprep_dir/", "文件夹", "PE去接头后fastq文件输出目录"],
        #     [r"./sickle_dir/", "文件夹", "质量剪切后fastq文件输出目录"]
        # ])
        super(MetagenomicQcModule, self).end()

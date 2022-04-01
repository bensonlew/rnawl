# -*- coding: utf-8 -*-
# __author__ = "wangzhaoyue"

import os
import re
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from collections import defaultdict


class FastpModule(Module):
    """
    Fastp质控
    author: wangzhaoyue
    last_modify: 2017.12.13
    """
    def __init__(self, work_id):
        super(FastpModule, self).__init__(work_id)
        options = [
            {"name": "sample_path", "type": "infile", "format": "datasplit.list_file"},
            # fastq路径list.txt文件，第一列路径，第二列样本名，第三列序列类型 l or r
            {"name": "qualified_quality_phred", "type": "string", "default": "20"},
            # -q,一个碱基合格的质量值,默认表示phred质量> = Q是合格的。
            {"name": "length_required", "type": "string", "default": "36"},   # -l,长度过滤参数，比此值短的读取将被丢弃
            {"name": "cut_by_quality5", "type": "string", "default": "20"},   # -5,根据前面(5 ')的质量，允许每个读切割，默认是禁用的
            {"name": "cut_by_quality3", "type": "string", "default": "3"},    # -3,根据后面(3 ')的质量，允许每个读切割，默认是禁用的
            {"name": "cut_mean_quality", "type": "string", "default": "20"},  # -M,在滑动窗口的基础上的平均质量低于切割质量将被切割，默认是Q20
            {"name": "n_base_limit", "type": "string", "default": "10"},      # -n,如果reads的碱基数大于该值，那么这个reads就被丢弃了
            {"name": "compression", "type": "string", "default": "6"},        # -z,gzip输出的压缩级别(1 ~ 9). 1是最快的，9是最小的
            {"name": "thread", "type": "string", "default": "8"},             # -w,线程数
            {"name": "cut_window_size", "type": "string"},                    # -W
            {"name": "adapter_sequence", "type": "string" },  # --adapter_sequence,the adapter for read1
            {"name": "adapter_sequence_r2", "type": "string"},  # --adapter_sequence_r2,the adapter for read2 (PE data only)
            # {"name": "adapter_sequence", "type": "string", "default": "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"},  # --adapter_sequence,the adapter for read1
            # {"name": "adapter_sequence_r2", "type": "string", "default": "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"},  # --adapter_sequence_r2,the adapter for read2 (PE data only)
        ]
        self.sample_info = defaultdict(list)
        self.tools = []
        self.end_times = 0
        self.add_option(options)

    def check_options(self):
        if not self.option("sample_path").is_set:
            raise OptionError("必须输入文库文件夹对应的路径信息")

    def run_fastp(self):
        for sample in self.sample_info:
            self.fastp = self.add_tool("datasplit_v2.fastp")
            opts = {
                "fq1": self.sample_info[sample][0],
                "fq2": self.sample_info[sample][1],
                "qualified_quality_phred": self.option("qualified_quality_phred"),
                "length_required": self.option("length_required"),
                "cut_mean_quality": self.option("cut_mean_quality"),
                "n_base_limit": self.option("n_base_limit"),
                "compression": self.option("compression"),
                "thread": self.option("thread"),
                "cut_window_size": self.option("cut_window_size"),
                "adapter_sequence": self.option("adapter_sequence"),
                "adapter_sequence_r2": self.option("adapter_sequence_r2")
            }
            if self.sample_lib_type[sample]:
                if re.search("微量", self.sample_lib_type[sample]):  # 真核RNA微量建库 接头read1和read2是一样的，都是CTGTCTCTTATACACATC
                    opts["adapter_sequence"] = "CTGTCTCTTATACACATC"
                    opts["adapter_sequence_r2"] = "CTGTCTCTTATACACATC"
            if len(self.sample_info[sample]) == 4:
                if self.sample_info[sample][2]:
                    opts["adapter_sequence"] = self.sample_info[sample][2]
                if self.sample_info[sample][3]:
                    opts["adapter_sequence_r2"] = self.sample_info[sample][3]
            if self.option("cut_by_quality5"):
                opts.update({"cut_by_quality5": self.option("cut_by_quality5")})
            if self.option("cut_by_quality3"):
                opts.update({"cut_by_quality3": self.option("cut_by_quality3")})
            self.fastp.set_options(opts)
            self.fastp.on("end", self.set_output, "fastp_{}".format(sample))  # modified by zengjing 2017.12.28(目的输出文件以样本命名)
            self.tools.append(self.fastp)
        for tool in self.tools:
            tool.run()

    def run(self):
        super(FastpModule, self).run()
        self.get_info()
        self.run_fastp()

    def set_output(self, event):
        obj = event["bind_object"]
        s = event["data"].split("fastp_")[1]
        json_dir = os.path.join(self.output_dir, "json_dir")
        fastq_dir = os.path.join(self.output_dir, "fastq")
        if not os.path.exists(json_dir):
            os.mkdir(json_dir)
        if not os.path.exists(fastq_dir):
            os.mkdir(fastq_dir)
        for f in os.listdir(obj.output_dir):
            m = re.match(r".+(.clean.1.fastq.*)", f)
            n = re.match(r".+(.clean.2.fastq.*)", f)
            if m:
                f2 = os.path.join(fastq_dir, s + m.group(1))
            elif n:
                f2 = os.path.join(fastq_dir, s + n.group(1))
            else:
                f2 = os.path.join(json_dir, s + ".json")
            if os.path.exists(f2):
                os.remove(f2)
            os.link(os.path.join(obj.output_dir, f), f2)
        self.end_times += 1
        if len(self.tools) == self.end_times:
            self.end()

    def end(self):
        super(FastpModule, self).end()

    def get_info(self):
        self.sample_lib_type = {}
        sample_size_null = {}
        with open(self.option("sample_path").prop["path"]) as fr:
            for line in fr:
                tmp = line.strip().split("\t")
                if len(tmp) < 3:
                    raise Exception("进行fastp的sample_path：{}文件必须大于三列，请检查！".format(self.option("sample_path").prop["path"]))
                size_null = self.check_fq_size(tmp[0])
                if size_null:
                    sample_size_null[tmp[1]] = size_null
                    continue
                if tmp[1] in sample_size_null.keys() and sample_size_null[tmp[1]]:
                        continue
                if tmp[1] in self.sample_info.keys():
                    if tmp[2] == "l":
                        self.sample_info[tmp[1]].insert(0, tmp[0])
                    else:
                        self.sample_info[tmp[1]].append(tmp[0])
                else:
                    self.sample_info[tmp[1]].append(tmp[0])
                if len(tmp) == 6:
                    size_null = self.check_fq_size(tmp[4])
                    if size_null:
                        sample_size_null[tmp[1]] = size_null
                        continue
                    size_null = self.check_fq_size(tmp[5])
                    if size_null:
                        sample_size_null[tmp[1]] = size_null
                        continue
                    if tmp[1] in sample_size_null.keys() and sample_size_null[tmp[1]]:
                        continue
                    self.sample_info[tmp[1]].append(tmp[4])
                    self.sample_info[tmp[1]].append(tmp[5])
                if len(tmp) > 3:
                    self.sample_lib_type[tmp[1]] = tmp[3]
            for key in self.sample_info.keys():
                if len(self.sample_info[key]) > 2:
                    raise Exception("需要质控的序列样本{}有重名，请改样本名或分开质控！".format(key))
                elif len(self.sample_info[key]) < 2:
                    raise Exception("样本{}对应的R1,R2序列不全,请核实！".format(key))

    def check_fq_size(self, path):
        size_null = False
        if os.path.isfile(path):
            lines_num = 0
            with open(path) as ft:
                for line_t in ft:
                    lines_num = lines_num + 1
                    if lines_num > 4:
                        break
            if lines_num < 4:
                self.logger.info("文件%s是空的，跳过这个样本的质控" % path)
                size_null = True
        return size_null

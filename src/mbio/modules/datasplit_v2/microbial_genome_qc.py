# -*- coding: utf-8 -*-
# __author__ = "zengjing"
# last_modify: 20181227

import os
import re
from biocluster.core.exceptions import OptionError
from biocluster.module import Module


class MicrobialGenomeQcModule(Module):
    """
    将样本去phix之后用fastp质控
    """
    def __init__(self, work_id):
        super(MicrobialGenomeQcModule, self).__init__(work_id)
        options = [
            {"name": "sample_info", "type": "infile", "format": "datasplit.micro_sample_info"},  # 样本信息表
            {"name": "flag", "type": "string", "default": "4"},  # 提取没有比对上的reads,此处固定取值4，软件默认0
            {"name": "readl", "type": "string"},  # 切除序列的阈值
            {"name": "qualified_quality_phred", "type": "string", "default": "20"},
            # -q,一个碱基合格的质量值,默认表示phred质量> = Q是合格的。
            {"name": "length_required", "type": "string", "default": "36"},  # -l,长度过滤参数，比此值短的读取将被丢弃
            {"name": "cut_by_quality5", "type": "string", "default": "20"},  # -5,根据前面(5 ")的质量，允许每个读切割，默认是禁用的
            {"name": "cut_by_quality3", "type": "string", "default": "3"},  # -3,根据后面(3 ")的质量，允许每个读切割，默认是禁用的
            {"name": "cut_mean_quality", "type": "string", "default": "20"},  # -M,在滑动窗口的基础上的平均质量低于切割质量将被切割，默认是Q20
            {"name": "n_base_limit", "type": "string", "default": "10"},  # -n,如果reads的碱基数大于该值，那么这个reads就被丢弃了
            {"name": "compression", "type": "string", "default": "6"},  # -z,gzip输出的压缩级别(1 ~ 9). 1是最快的，9是最小的
            {"name": "thread", "type": "string", "default": "8"},  # -w,线程数
            {"name": "adapter_sequence", "type": "string", "default": "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"},  # --adapter_sequence,the adapter for read1
            {"name": "adapter_sequence_r2", "type": "string", "default": "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"},  # --adapter_sequence_r2,the adapter for read2 (PE data only)
        ]
        self.sample_info = {}
        self.all_samples = []
        self.end_times = 0
        self.add_option(options)

    def check_options(self):
        if not self.option("sample_info").is_set:
            raise OptionError("必须输入样本信息文件，包括样本，文库类型、插入片段长度三列")

    def get_info(self):
        """
        获得样本对应的路径信息，以及样本的
        """
        with open(self.option("sample_info").prop["path"])as f:
            lines = f.readlines()
            for line in lines[1:]:
                tmp = line.strip().split("\t")
                name = tmp[1] + "--" + tmp[0]
                self.sample_info[name] = {"lib_size": tmp[2], "path": tmp[4].split(";")}
                if len(tmp) == 7:
                    self.sample_info[name]["adapter1"] = tmp[5]
                    self.sample_info[name]["adapter1"] = tmp[6]
                self.all_samples.append(name)
        for s in self.all_samples:
            if len(self.sample_info[s]["path"]) != 2:
                raise Exception("样本{}对应的R1,R2序列有误,请检查！".format(s))

    def run_phix_filter(self):
        for s in self.all_samples:
            options = {
                "fq1": self.sample_info[s]["path"][0],
                "fq2": self.sample_info[s]["path"][1],
                "sample_name": s,
                "insert_size": self.sample_info[s]["lib_size"],
                "flag": self.option("flag"),
            }
            if self.option("readl"):
                options["readl"] = self.option("readl")
            self.phix_filter = self.add_tool("datasplit_v2.phix_filter")
            self.phix_filter.set_options(options)
            self.phix_filter.on("end", self.run_fastp, s)
            self.phix_filter.run()

    def run_fastp(self, event):
        obj = event["bind_object"]
        options = {
            "fq1": obj.output_dir + "/" + event["data"] + ".filtphix.1.result.fq",
            "fq2": obj.output_dir + "/" + event["data"] + ".filtphix.2.result.fq",
            "qualified_quality_phred": self.option("qualified_quality_phred"),
            "length_required": self.option("length_required"),
            "cut_mean_quality": self.option("cut_mean_quality"),
            "n_base_limit": self.option("n_base_limit"),
            "compression": self.option("compression"),
            "thread": self.option("thread"),
            "adapter_sequence": self.option("adapter_sequence"),
            "adapter_sequence_r2": self.option("adapter_sequence_r2")
        }
        if self.option("cut_by_quality5"):
            options["cut_by_quality5"] = self.option("cut_by_quality5")
        if self.option("cut_by_quality3"):
            options["cut_by_quality3"] = self.option("cut_by_quality3")
        if "adapter1" in self.sample_info[event["data"]].keys():
            if "adapter1" in self.sample_info[event["data"]].keys() and self.sample_info[event["data"]]["adapter1"]:
                options["adapter_sequence"] = self.sample_info[event["data"]]["adapter1"]
            if "adapter2" in self.sample_info[event["data"]].keys() and self.sample_info[event["data"]]["adapter2"]:
                options["adapter_sequence_r2"] = self.sample_info[event["data"]]["adapter2"]
        self.fastp = self.add_tool("datasplit_v2.fastp")
        self.fastp.set_options(options)
        self.fastp.on("end", self.set_output, event["data"])
        self.fastp.run()

    def run(self):
        super(MicrobialGenomeQcModule, self).run()
        self.get_info()
        self.run_phix_filter()

    def set_output(self, event):
        self.logger.info("设置结果目录")
        obj = event["bind_object"]
        json_dir = os.path.join(self.output_dir, "json_dir")
        fastq_dir = os.path.join(self.output_dir, "fastq")
        if not os.path.exists(json_dir):
            os.mkdir(json_dir)
        if not os.path.exists(fastq_dir):
            os.mkdir(fastq_dir)
        for f in os.listdir(obj.output_dir):
            f1 = os.path.join(obj.output_dir, f)
            m = re.match(r".+(.clean.1.fastq.*)", f)
            n = re.match(r".+(.clean.2.fastq.*)", f)
            if m:
                f2 = os.path.join(fastq_dir, event["data"] + m.group(1))
            elif n:
                f2 = os.path.join(fastq_dir, event["data"] + n.group(1))
            else:
                f2 = os.path.join(json_dir, event["data"] + ".json")
            if os.path.exists(f2):
                os.remove(f2)
            os.link(f1, f2)
        self.end_times += 1
        if self.end_times == len(self.all_samples):
            self.end()

    def end(self):
        super(MicrobialGenomeQcModule, self).end()

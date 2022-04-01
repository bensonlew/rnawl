# -*- coding: utf-8 -*-

import os
import re
from collections import defaultdict
from biocluster.module import Module
from biocluster.core.exceptions import OptionError


class FastpModule(Module):
    """
    Fastp质控
    """
    def __init__(self, work_id):
        super(FastpModule, self).__init__(work_id)
        options = [
            {"name": "sample_list", "type": "infile", "format": "noref_wgs.list_file", "required": True},
            # 第一列分析名称，第二列批次样本名，第三列左端序列，第四列右端序列
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
        self.end_times = 0
        self.add_option(options)

    def check_options(self):
        pass

    def get_info(self):
        self.sample_list = []
        self.sample_info = {}
        with open(self.option("sample_list").prop["path"])as fr:
            for line in fr:
                tmp = line.strip().split("\t")
                if len(tmp) != 4:
                    self.set_error("进行fastp的sample_list：%s文件必须是四列，请检查！", variables=(self.option("sample_list").prop["path"]), code="25500609")
                if tmp[1] in self.sample_list:
                    self.set_error("sample_list:%s里的样本：%s重复，请检查", variables=(self.option("sample_list").prop["path"], tmp[1]), code="25500610")
                else:
                    self.sample_info[tmp[1]] = {"fq": [], "analysis_name": ""}
                    self.sample_info[tmp[1]]["analysis_name"] = tmp[0]
                    self.sample_info[tmp[1]]["fq"].append(tmp[2])
                    self.sample_info[tmp[1]]["fq"].append(tmp[3])
                    self.sample_list.append(tmp[1])

    def run_fastp(self):
        self.tools = []
        for sample in self.sample_list:
            self.fastp = self.add_tool("wgs.fastp")
            options = {
                "fq1": self.sample_info[sample]["fq"][0],
                "fq2": self.sample_info[sample]["fq"][1],
                "qualified_quality_phred": self.option("qualified_quality_phred"),
                "length_required": self.option("length_required"),
                "cut_mean_quality": self.option("cut_mean_quality"),
                "n_base_limit": self.option("n_base_limit"),
                "compression": self.option("compression"),
                "thread": self.option("thread"),
                "adapter_sequence": self.option("adapter_sequence"),
                "adapter_sequence_r2": self.option("adapter_sequence_r2"),
                "cut_by_quality3": self.option("cut_by_quality3"),
                "cut_by_quality5": self.option("cut_by_quality5"),
                "sample_name": sample
            }
            self.fastp.set_options(options)
            self.fastp.on("end", self.set_output, "fastp_{}".format(sample))
            self.tools.append(self.fastp)
        if self.tools:
            for tool in self.tools:
                tool.run()
        else:
            self.set_error("tool 队列是空的，无法进行后续的分析，流程终止", code="25500611")

    def set_output(self, event):
        self.logger.info("设置结果目录")
        obj = event["bind_object"]
        self.end_times += 1
        m = re.match(r"fastp_(.+)", event["data"])
        s = m.group(1)
        json_dir = os.path.join(self.output_dir, "json_dir")
        fastq_dir = os.path.join(self.output_dir, "fastq")
        if not os.path.exists(json_dir):
            os.mkdir(json_dir)
        if not os.path.exists(fastq_dir):
            os.mkdir(fastq_dir)
        json_path = obj.work_dir + "/" + s + ".json"
        fq1_path = obj.output_dir + "/" + s + ".clean.1.fastq.gz"
        fq2_path = obj.output_dir + "/" + s + ".clean.2.fastq.gz"
        fq1_path_ = os.path.join(fastq_dir, s + ".clean.1.fastq.gz")
        fq2_path_ = os.path.join(fastq_dir, s + ".clean.2.fastq.gz")
        self.link_file(json_path, os.path.join(json_dir, s + ".json"))
        self.link_file(fq1_path, fq1_path_)
        self.link_file(fq2_path, fq2_path_)
        if len(self.tools) == self.end_times:
            with open(os.path.join(self.output_dir, "fastq.list"), "w") as w:
                for s in self.sample_list:
                    fq1_path = os.path.join(fastq_dir, s + ".clean.1.fastq.gz")
                    fq2_path = os.path.join(fastq_dir, s + ".clean.2.fastq.gz")
                    w.write(self.sample_info[s]["analysis_name"] + "\t" + s + "\t" + fq1_path + "\t" + fq2_path + "\n")
            self.end()

    def link_file(self, old, new):
        if not os.path.exists(old):
            self.set_error("需要进行链接的原始文件：%s 不存在，请检查" , variables=( old), code="25500612")
        if os.path.exists(new):
            os.remove(new)
        os.link(old, new)

    def run(self):
        super(FastpModule, self).run()
        self.get_info()
        self.run_fastp()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        super(FastpModule, self).end()

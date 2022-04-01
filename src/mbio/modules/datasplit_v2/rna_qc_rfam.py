# -*- coding: utf-8 -*-
# __author__ = "zengjing"

import os
import re
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from collections import defaultdict


class RnaQcRfamModule(Module):
    """
    Fastp质控 + blast Rfam统计，用于lncrna和prokaryotic_rna
    """
    def __init__(self, work_id):
        super(RnaQcRfamModule, self).__init__(work_id)
        options = [
            {"name": "sample_path", "type": "infile", "format": "datasplit.list_file"},
            # fastq路径list.txt文件，第一列路径，第二列样本名，第三列序列类型 l or r
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
            {"name": "project_type", "type": "string", "default": "lncrna"},
            {"name": "database", "type": "string", "default": "rfam"},  # 数据库，nt或者rfam
        ]
        self.sample_info = defaultdict(list)
        self.end_times = 0
        self.add_option(options)

    def check_options(self):
        if not self.option("sample_path").is_set:
            raise OptionError("必须输入文库文件夹对应的路径信息")
        if self.option("project_type") not in ["lncrna", "prokaryotic_rna"]:
            raise OptionError("项目类型只能是lncrna和prokaryotic_rna")
        if self.option("database") not in ["nt", "rfam"]:
            raise OptionError("项目类型只能是nt和rfam")

    def get_info(self):
        with open(self.option("sample_path").prop["path"])as fr:
            for line in fr:
                tmp = line.strip().split("\t")
                if len(tmp) < 3:
                    raise Exception("进行fastp的sample_path：{}文件必须是三列，请检查！".format(self.option("sample_path").prop["path"]))
                if tmp[1] in self.sample_info.keys():
                    if tmp[2] == "l":
                        self.sample_info[tmp[1]].insert(0, tmp[0])
                    else:
                        self.sample_info[tmp[1]].append(tmp[0])
                else:
                    self.sample_info[tmp[1]].append(tmp[0])
                if len(tmp) == 6:
                    self.sample_info[tmp[1]].append(tmp[4])
                    self.sample_info[tmp[1]].append(tmp[5])
            for key in self.sample_info.keys():
                if len(self.sample_info[key]) > 2:
                    raise Exception("需要质控的序列样本{}有重名，请改样本名或分开质控！".format(key))
                elif len(self.sample_info[key]) < 2:
                    raise Exception("样本{}对应的R1,R2序列不全,请核实！".format(key))

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
                "adapter_sequence": self.option("adapter_sequence"),
                "adapter_sequence_r2": self.option("adapter_sequence_r2")
            }
            if self.option("cut_by_quality5"):
                opts.update({"cut_by_quality5": self.option("cut_by_quality5")})
            if self.option("cut_by_quality3"):
                opts.update({"cut_by_quality3": self.option("cut_by_quality3")})
            if len(self.sample_info[sample]) == 4:
                if self.sample_info[sample][2]:
                    opts["adapter_sequence"] = self.sample_info[sample][2]
                if self.sample_info[sample][3]:
                    opts["adapter_sequence_r2"] = self.sample_info[sample][3]
            self.fastp.set_options(opts)
            self.fastp.on("end", self.set_output, "fastp:{}".format(sample))
            # self.fastp.on("end", self.run_blast_rfam, sample)
            self.fastp.run()

    def run_blast_rfam(self, event):
        """
        随机抽取10000条序列进行Rfam比对，统计核糖体占比
        """
        obj = event["bind_object"]
        for f in os.listdir(obj.output_dir):
            m = re.match(r".+(.clean.1.fastq.*)", f)
            if m:
                fastq = os.path.join(obj.output_dir, f)
        options = {
            "fastq": fastq,
            "database": self.option("database"),
            "project_type": self.option("project_type")
        }
        self.blast_rfam = self.add_module("datasplit.single_fastq_blast")
        self.blast_rfam.set_options(options)
        self.blast_rfam.on("end", self.set_output, "rfam:{}".format(event["data"]))
        self.blast_rfam.run()

    def run(self):
        super(RnaQcRfamModule, self).run()
        self.get_info()
        self.run_fastp()

    def link(self, old, new):
        if os.path.exists(new):
            os.remove(new)
        os.link(old, new)

    def set_output(self, event):
        obj = event["bind_object"]
        rfam_dir = os.path.join(self.output_dir, "rfam_dir")
        json_dir = os.path.join(self.output_dir, "json_dir")
        fastq_dir = os.path.join(self.output_dir, "fastq")
        if not os.path.exists(rfam_dir):
            os.mkdir(rfam_dir)
        if not os.path.exists(fastq_dir):
            os.mkdir(fastq_dir)
        if not os.path.exists(json_dir):
            os.mkdir(json_dir)
        if event["data"].startswith("fastp"):
            s = event["data"].split("fastp:")[1]
            for f in os.listdir(obj.output_dir):
                m = re.match(r".+(.clean.1.fastq.*)", f)
                n = re.match(r".+(.clean.2.fastq.*)", f)
                if m:
                    f2 = os.path.join(fastq_dir, s + m.group(1))
                elif n:
                    f2 = os.path.join(fastq_dir, s + n.group(1))
                else:
                    f2 = os.path.join(json_dir, s + ".json")
                self.link(os.path.join(obj.output_dir, f), f2)
        elif event["data"].startswith("rfam"):
            self.link(os.path.join(obj.output_dir, "rfam_summary.xls"), os.path.join(rfam_dir, event["data"].split("rfam:")[1] + ".rfam_summary.xls"))
        self.end_times += 1
        # if 2 * len(self.sample_info) == self.end_times:
        if len(self.sample_info) == self.end_times:
            self.end()

    def end(self):
        super(RnaQcRfamModule, self).end()

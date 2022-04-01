# -*- coding: utf-8 -*-
# __author__ = "zengjing"
# modified 20190305

import os
import re
from biocluster.module import Module
from biocluster.core.exceptions import OptionError


class QcStatModule(Module):
    """
    WGS项目fastq序列用fastp软件质控，用ngsqc软件统计质控前后的序列信息
    """
    def __init__(self, work_id):
        super(QcStatModule, self).__init__(work_id)
        options = [
            {"name": "fastq_dir", "type": "infile", "format": "sequence.fastq_dir"},  # 输入的fq文件夹
            {"name": "task_id", "type": "string"},  # 任务id,从sg_specimen_other表和sg_task表得到fastq_list
            {"name": "qualified_quality_phred", "type": "string", "default": "20"},  # -q,一个碱基合格的质量值,默认表示phred质量> = Q是合格的。
            {"name": "length_required", "type": "string", "default": "36"},  # -l,长度过滤参数，比此值短的读取将被丢弃
            {"name": "cut_by_quality5", "type": "string", "default": "20"},  # -5,根据前面(5 ')的质量，允许每个读切割，默认是禁用的
            {"name": "cut_by_quality3", "type": "string", "default": "3"},  # -3,根据后面(3 ')的质量，允许每个读切割，默认是禁用的
            {"name": "cut_mean_quality", "type": "string", "default": "20"},  # -M,在滑动窗口的基础上的平均质量低于切割质量将被切割，默认是Q20
            {"name": "n_base_limit", "type": "string", "default": "10"},  # -n,如果reads的碱基数大于该值，那么这个reads就被丢弃了
            {"name": "compression", "type": "string", "default": "6"},  # -z,gzip输出的压缩级别(1 ~ 9). 1是最快的，9是最小的
            {"name": "thread", "type": "string", "default": "8"},  # -w,线程数
            {"name": "adapter_sequence", "type": "string", "default": "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"},  # --adapter_sequence,the adapter for read1
            {"name": "adapter_sequence_r2", "type": "string", "default": "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"},  # --adapter_sequence_r2,the adapter for read2 (PE data only)
            {"name": "project_type", "type": "string", "default": "dna_wgs_v2"}
        ]
        self.add_option(options)
        self.clean_data = {}
        self.times = 0

    def check_options(self):
        if not self.option("task_id"):
            raise OptionError("请设置task_id", code="24501101")

    def get_specimen_info(self):
        """
        从mongo库的sg_task、sg_specimen_other表中得到样本及对应的fastq
        """
        qc_fastq = self.api.api("wgs.qc_fastq")
        if self.option("project_type") == "dna_wgs":
            pass
        else:
            qc_fastq._project_type = self.option("project_type")
        if self.option("fastq_dir").is_set:
            self.specimen_info, self.rawdata_info = qc_fastq.get_specimen_fastq(self.option("task_id"), self.option("fastq_dir").prop["path"])
        else:
            self.specimen_info, self.rawdata_info = qc_fastq.get_specimen_fastq(self.option("task_id"))

    def run_fastq_qc(self):
        """
        质控
        """
        for s in self.rawdata_info.keys():
            options = {
                "fq1": self.rawdata_info[s]["l"],
                "fq2": self.rawdata_info[s]["r"],
                "qualified_quality_phred": "20",
                "length_required": "36",
                "cut_by_quality5": "20",
                "cut_by_quality3": "3",
                "cut_mean_quality": "20",
                "n_base_limit": "10",
                "thread": "8",
                "compression": "6",
                "sample_name": s,
                "adapter_sequence": self.option("adapter_sequence"),
                "adapter_sequence_r2": self.option("adapter_sequence_r2")
            }
            self.fastq_qc = self.add_tool("wgs.fastp")
            self.fastq_qc.set_options(options)
            self.fastq_qc.on("end", self.set_output, s)
            self.fastq_qc.run()

    def link_file(self, old, new):
        if os.path.exists(new):
            os.remove(new)
        if os.path.exists(old):
            os.link(old, new)
        else:
            self.set_error("文件：%s不存在，请检查", variables=(old), code="24501107")

    def set_output(self, event):
        """
        得到质控后的clean.fastq文件的字典集合
        """
        obj = event["bind_object"]
        sample_name = event["data"]
        if sample_name not in self.clean_data.keys():
            self.clean_data[sample_name] = {}
        qc_dir = os.path.join(self.output_dir, "clean_data")
        json_dir = os.path.join(self.output_dir, "qc_stat")
        if not os.path.exists(qc_dir):
            os.mkdir(qc_dir)
        if not os.path.exists(json_dir):
            os.mkdir(json_dir)
        for f in os.listdir(obj.output_dir):
            self.link_file(obj.output_dir + "/" + f, qc_dir + "/" + f)
            if re.search(r".*clean.1.fastq.gz", f):
                self.clean_data[sample_name]["l"] = qc_dir + "/" + f
            if re.search(r".*clean.2.fastq.gz", f):
                self.clean_data[sample_name]["r"] = qc_dir + "/" + f
        for f in os.listdir(obj.work_dir):
            if f.endswith(".json"):
                self.link_file(obj.work_dir + "/" + f, json_dir + "/" + f)
        if len(self.clean_data.keys()) == len(self.rawdata_info.keys()):
            with open(os.path.join(qc_dir, "fastq.list"), "w") as w:
                for s in self.clean_data.keys():
                    w.write(s + "\t" + self.clean_data[s]["l"] + "\t" + self.clean_data[s]["r"] + "\t" + self.rawdata_info[s]["lib"] + "\n")
        self.times += 1
        if self.times == len(self.rawdata_info.keys()):
            self.end()

    def run(self):
        super(QcStatModule, self).run()
        self.get_specimen_info()
        self.logger.info(self.rawdata_info)
        if len(self.rawdata_info.keys()) == 0:
            self.set_error("没有找到对应的样本信息，请检查", code="24501109")
        self.run_fastq_qc()

    def end(self):
        super(QcStatModule, self).end()

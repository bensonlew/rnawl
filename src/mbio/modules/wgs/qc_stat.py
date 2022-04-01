# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# modified 2018.04.03

from biocluster.core.exceptions import OptionError
from biocluster.module import Module
import re
import os
import json
import time


class QcStatModule(Module):
    """
    WGS项目fastq序列用fastp软件质控，用ngsqc软件统计质控前后的序列信息
    """
    def __init__(self, work_id):
        super(QcStatModule, self).__init__(work_id)
        options = [
            {'name': 'fastq_dir', 'type': 'infile', 'format': 'sequence.fastq_dir'},  # 输入的fq文件夹
            # {"name": "fastq_list", "type": "infile", "format": "wgs.fastq_list"},
            # fastq路径list.txt文件，第一列样本名称，第二列左端序列，第三列右端序列，第四列文库类型
            {'name': 'task_id', 'type': 'string'},  # 任务id,从sg_specimen_other表和sg_task表得到fastq_list
            {'name': 'qualified_quality_phred', 'type': "string", "default": "20"},
            # -q,一个碱基合格的质量值,默认表示phred质量> = Q是合格的。
            {'name': 'length_required', "type": "string", "default": "36"},  # -l,长度过滤参数，比此值短的读取将被丢弃
            {'name': "cut_by_quality5", "type": "string", "default": "20"},  # -5,根据前面(5 ')的质量，允许每个读切割，默认是禁用的
            {'name': "cut_by_quality3", "type": "string", "default": "3"},  # -3,根据后面(3 ')的质量，允许每个读切割，默认是禁用的
            {'name': 'cut_mean_quality', "type": "string", "default": "20"},  # -M,在滑动窗口的基础上的平均质量低于切割质量将被切割，默认是Q20
            {'name': 'n_base_limit', "type": "string", "default": "10"},  # -n,如果reads的碱基数大于该值，那么这个reads就被丢弃了
            {'name': 'compression', "type": "string", "default": "6"},  # -z,gzip输出的压缩级别(1 ~ 9). 1是最快的，9是最小的
            {'name': 'thread', "type": "string", "default": "8"},  # -w,线程数
            {'name': "project_type", "type": "string", "default": "dna_wgs"}
        ]
        self.add_option(options)
        self.clean_data = {}
        self.times, self.raw_times = 0, 0

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
            qc_fastq.set_project_type_as_workflow()  # 设置api与workflow一致，add by hd@20180927
        if self.option("fastq_dir").is_set:
            self.specimen_info, self.rawdata_info = qc_fastq.get_specimen_fastq(self.option("task_id"), self.option("fastq_dir").prop["path"])
        else:
            self.specimen_info, self.rawdata_info = qc_fastq.get_specimen_fastq(self.option("task_id"))

    def run_raw_data_stat(self):
        """
        原始序列数据统计
        """
        for s in self.rawdata_info.keys():
            options = {
                "fastq_l": self.rawdata_info[s]["l"],
                "fastq_r": self.rawdata_info[s]["r"],
                "sample_name": s
            }
            self.rawdata_stat = self.add_tool("wgs.ngsqc_stat")
            self.rawdata_stat.set_options(options)
            self.rawdata_stat.on("end", self.set_output, "raw_{}".format(s))
            self.rawdata_stat.run()

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
                "sample_name": s
            }
            self.fastq_qc = self.add_tool("wgs.fastp")
            self.fastq_qc.set_options(options)
            self.fastq_qc.on("end", self.set_output, "fastq_qc:{}".format(s))
            self.fastq_qc.run()

    def run_clean_data_stat(self):
        """
        质控后序列数据统计
        """
        stat = []
        for s in self.clean_data.keys():
            options = {
                "fastq_l": self.clean_data[s]["l"],
                "fastq_r": self.clean_data[s]["r"],
                "sample_name": s
            }
            self.clean_data_stat = self.add_tool("wgs.ngsqc_stat")
            self.clean_data_stat.set_options(options)
            self.clean_data_stat.on("end", self.set_output, "clean_{}".format(s))
            self.clean_data_stat.run()

    def link_file(self, old, new):
        if os.path.exists(new):
            os.remove(new)
        if os.path.exists(old):
            os.link(old, new)
        else:
            self.set_error("文件：%s不存在，请检查", variables=(old), code="24501107")
        if os.path.exists(new):
            return True
        else:
            for i in range(5):
                time.sleep(1)
                self.link_file(old, new)

    def set_output(self, event):
        """
        得到质控后的clean.fastq文件的字典集合
        """
        obj = event["bind_object"]
        m = re.match(r"fastq_qc:(.+)", event["data"])
        n = re.match(r"raw_(.+)", event["data"])
        l = re.match(r"clean_(.+)", event["data"])
        if m:
            sample_name = m.group(1)
            clean_fastq = self.output_dir + "/clean_data"
            if not os.path.exists(clean_fastq):
                os.mkdir(clean_fastq)
            if sample_name not in self.clean_data.keys():
                self.clean_data[sample_name] = {}
            for f in os.listdir(obj.output_dir):
                self.link_file(obj.output_dir + "/" + f, clean_fastq + "/" + f)
                if re.search(r".*clean.1.fastq.gz", f):
                    self.clean_data[sample_name]["l"] = clean_fastq + "/" + f
                if re.search(r".*clean.2.fastq.gz", f):
                    self.clean_data[sample_name]["r"] = clean_fastq + "/" + f
            if len(self.clean_data.keys()) == len(self.rawdata_info.keys()):
                with open(os.path.join(clean_fastq, "fastq.list"), "w") as w:
                    for s in self.clean_data.keys():
                        w.write(s + "\t" + self.clean_data[s]["l"] + "\t" + self.clean_data[s]["r"] + "\t" + self.rawdata_info[s]["lib"] + "\n")
                self.run_clean_data_stat()
        if n:
            raw_stat = self.output_dir + "/rawdata_qc"
            raw_atgc = raw_stat + "/atgc"
            raw_qual = raw_stat + "/qual"
            sample_name = n.group(1)
            if not os.path.exists(raw_stat):
                os.mkdir(raw_stat)
            if not os.path.exists(raw_atgc):
                os.mkdir(raw_atgc)
            if not os.path.exists(raw_qual):
                os.mkdir(raw_qual)
            for f in os.listdir(obj.output_dir):
                if re.search(r".*.atgc", f):
                    self.link_file(obj.output_dir + "/" + f, raw_atgc + "/" + f + ".xls")
                if re.search(r".*.qual", f):
                    self.link_file(obj.output_dir + "/" + f, raw_qual + "/" + f + ".xls")
                if re.search(r".*.stat", f):
                    raw_qc = os.path.join(self.work_dir, "raw_qc_stat")
                    if not os.path.exists(raw_qc):
                        os.mkdir(raw_qc)
                    self.link_file(obj.output_dir + "/" + f, raw_qc + "/" + f)
            self.raw_times += 1
            if self.raw_times == len(self.rawdata_info.keys()):
                header = "#sampleID\tRaw Reads\tRaw Bases(bp)\tRaw GC(%)\tRaw Q30(%)\n"
                self.get_all_specimen_stat(header, os.path.join(self.work_dir, "raw_qc_stat"), raw_stat + "/qc.xls")
        if l:
            clean_stat = self.output_dir + "/cleandata_qc"
            clean_atgc = clean_stat + "/atgc"
            clean_qual = clean_stat + "/qual"
            sample_name = l.group(1)
            if not os.path.exists(clean_stat):
                os.mkdir(clean_stat)
            if not os.path.exists(clean_atgc):
                os.mkdir(clean_atgc)
            if not os.path.exists(clean_qual):
                os.mkdir(clean_qual)
            for f in os.listdir(obj.output_dir):
                if re.search(r".*.atgc", f):
                    self.link_file(obj.output_dir + "/" + f, clean_atgc + "/" + f + ".xls")
                if re.search(r".*.qual", f):
                    self.link_file(obj.output_dir + "/" + f, clean_qual + "/" + f + ".xls")
                if re.search(r".*.stat", f):
                    clean_qc = os.path.join(self.work_dir, "clean_qc_stat")
                    if not os.path.exists(clean_qc):
                        os.mkdir(clean_qc)
                    self.link_file(obj.output_dir + "/" + f, clean_qc + "/" + f)
            self.times += 1
            if self.times == len(self.rawdata_info.keys()):
                header = "#sampleID\tClean Reads\tClean Bases(bp)\tClean GC(%)\tClean Q30(%)\n"
                self.get_all_specimen_stat(header, os.path.join(self.work_dir, "clean_qc_stat"), clean_stat + "/qc.xls")
                self.end()

    def get_all_specimen_stat(self, header, stat_dir, qc_stat):
        """
        将多个样本的stat写入一个文件
        """
        # time.sleep(10)
        with open(qc_stat, "a") as w:
            w.write(header)
            for analysis_name in self.specimen_info.keys():
                reads, base, gc, q30 = 0, 0, 0, 0
                for name in self.specimen_info[analysis_name]:
                    sample_stat = os.path.join(stat_dir, name + ".stat")
                    if not os.path.exists(sample_stat):
                        self.set_error("没有找到样本%s对应的stat文件%s", variables=(name, sample_stat), code="24501108")
                    with open(sample_stat, "r") as f:
                        lines = f.readlines()
                        item = lines[-1].strip().split("\t")
                        reads += int(item[1])
                        base += int(item[2])
                        gc += float(item[3])
                        q30 += float(item[4])
                len_ = len(self.specimen_info[analysis_name])
                gc = round(float(gc) / len_, 2)
                q30 = round(float(q30) / len_, 2)
                w.write(analysis_name + "\t" + str(reads/2) + "\t" + str(base) + "\t" + str(gc) + "\t" +
                        str(q30) + "\n")

    def run(self):
        super(QcStatModule, self).run()
        raw_qc = self.output_dir + "/rawdata_qc/qc.xls"
        clean_qc = self.output_dir + "/cleandata_qc/qc.xls"
        if os.path.exists(raw_qc):
            os.remove(raw_qc)
        if os.path.exists(clean_qc):
            os.remove(clean_qc)
        self.get_specimen_info()
        self.logger.info(self.rawdata_info)
        if len(self.rawdata_info.keys()) == 0:
            self.set_error("没有找到对应的样本信息，请检查", code="24501109")
        self.run_raw_data_stat()
        self.run_fastq_qc()

    def end(self):
        super(QcStatModule, self).end()

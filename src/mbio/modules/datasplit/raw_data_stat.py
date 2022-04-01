# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# last_modify: 20171214

import os
import re
import shutil
from biocluster.module import Module
from biocluster.core.exceptions import OptionError


class RawDataStatModule(Module):
    """
    原始数据统计，统计冗余度，fastq_dup.py(输入为文件夹时文件夹内需要有list.txt)
    统计碱基质量信息及Q20、Q30，fastx_quality_stats、q20q30_stat.py(输入为文件夹时文件夹内需要有list.txt)
    统计fastq序列基本信息，FastqStat.jar(输入为文件夹时文件夹内需要有list.txt)
    """
    def __init__(self, work_id):
        super(RawDataStatModule, self).__init__(work_id)
        options = [
            {"name": "list_file", "type": "infile", "format": "datasplit.list_file"},  # 要进行质控的文件的list,第一列文件路径，二列文库，三列r/l
        ]
        self.add_option(options)
        self.step.add_steps("dup", "stat", "fastx")
        self.samples = {}
        self.start_times = 0
        self.end_times = 0

    def set_step(self, event):
        if "start" in event["data"].keys():
            event["data"]["start"].start()
        if "end" in event["data"].keys():
            event["data"]["end"].finish()
        self.step.update()

    def check_options(self):
        if not self.option("list_file").is_set:
            raise OptionError("必须设置进行质控的文件list")
        self.samples = self.option("list_file").prop["samples"]

    def run_fastq_dup(self):
        for s in self.samples.keys():
            if len(self.samples[s]) == 1:
                options = {
                    "fastq_s": self.samples[s][0],
                    "fq_type": "SE"
                }
            else:
                options = {
                    "fastq_l": self.samples[s][0],
                    "fastq_r": self.samples[s][1],
                    "fq_type": "PE"
                }
            self.fastq_dup = self.add_tool("datasplit.fastq_dup")
            self.fastq_dup.set_options(options)
            self.fastq_dup.on("start", self.set_step, {"start": self.step.dup})
            self.fastq_dup.on("end", self.set_step, {"end": self.step.dup})
            self.fastq_dup.on('end', self.set_output, 'fastq_dup_{}'.format(s))
            self.fastq_dup.run()
            self.start_times += 1

    def run_fastq_stat(self):
        options = {
            "list_file": self.option("list_file")
        }
        self.fastq_stat = self.add_tool("datasplit.fastq_stat")
        self.fastq_stat.set_options(options)
        self.fastq_stat.on("start", self.set_step, {"start": self.step.stat})
        self.fastq_stat.on("end", self.set_step, {"end": self.step.stat})
        self.fastq_stat.on('end', self.set_output, 'fastq_stat')
        self.fastq_stat.run()
        self.start_times += 1

    def run_fastx(self):
        for s in self.samples.keys():
            if len(self.samples[s]) == 2:
                options = {
                    "fastq": self.samples[s][1]
                }
                self.fastx = self.add_tool("datasplit.fastx_v2")
                self.fastx.set_options(options)
                self.fastx.on("start", self.set_step, {"start": self.step.fastx})
                self.fastx.on("end", self.set_step, {"end": self.step.fastx})
                self.fastx.on('end', self.set_output, 'fastx_{}_r'.format(s))
                self.fastx.run()
                self.start_times += 1
                options = {
                    "fastq": self.samples[s][0]
                }
                self.fastx = self.add_tool("datasplit.fastx_v2")
                self.fastx.set_options(options)
                self.fastx.on("start", self.set_step, {"start": self.step.fastx})
                self.fastx.on("end", self.set_step, {"end": self.step.fastx})
                self.fastx.on('end', self.set_output, 'fastx_{}_l'.format(s))
                self.fastx.run()
                self.start_times += 1
            else:
                options = {
                    "fastq": self.samples[s][0]
                }
                self.fastx = self.add_tool("datasplit.fastx_v2")
                self.fastx.set_options(options)
                self.fastx.on("start", self.set_step, {"start": self.step.fastx})
                self.fastx.on("end", self.set_step, {"end": self.step.fastx})
                self.fastx.on('end', self.set_output, 'fastx_{}'.format(s))
                self.fastx.run()
                self.start_times += 1

    def set_output(self, event):
        self.logger.info("开始set output")
        obj = event["bind_object"]
        self.end_times += 1
        m1 = re.match(r"fastq_dup_(.+)", event["data"])
        m2 = re.match(r"fastx_(.+)", event["data"])
        if m1:
            new_dir = os.path.join(self.output_dir, "fastq_dup")
            if not os.path.exists(new_dir):
                os.mkdir(new_dir)
            for f in os.listdir(obj.output_dir):
                if os.path.exists(new_dir + "/" + m1.group(1) + "_dup.xls"):
                    os.remove(new_dir + "/" + m1.group(1) + "_dup.xls")
                os.link(obj.output_dir + "/" + f, new_dir + "/" + m1.group(1) + "_dup.xls")
        if m2:
            new_dir = os.path.join(self.output_dir, "fastx")
            if not os.path.exists(new_dir):
                os.mkdir(new_dir)
            for f in os.listdir(obj.output_dir):
                n1 = re.match(r".+fastxstat$", f)
                n2 = re.match(r".+q20q30$", f)
                if n1:
                    if os.path.exists(new_dir + "/" + m2.group(1) + "_fastxstat"):
                        os.remove(new_dir + "/" + m2.group(1) + "_fastxstat")
                    os.link(obj.output_dir + "/" + f, new_dir + "/" + m2.group(1) + "_fastxstat")
                if n2:
                    if os.path.exists(new_dir + "/" + m2.group(1) + "_q20q30"):
                        os.remove(new_dir + "/" + m2.group(1) + "_q20q30")
                    os.link(obj.output_dir + "/" + f, new_dir + "/" + m2.group(1) + "_q20q30")
        if event["data"] == "fastq_stat":
            new_dir = os.path.join(self.output_dir, "fastq_stat")
            if os.path.exists(new_dir):
                shutil.rmtree(new_dir)
            os.mkdir(new_dir)
            for f in os.listdir(obj.output_dir):
                if os.path.exists(new_dir + "/" + f):
                    os.remove(new_dir + "/" + f)
                os.link(obj.output_dir + "/" + f, new_dir + "/" + f)
        if self.start_times == self.end_times:
            self.logger.info("结束set output")
            self.end()

    def run(self):
        super(RawDataStatModule, self).run()
        self.run_fastq_dup()
        self.run_fastq_stat()
        self.run_fastx()

    def end(self):
        super(RawDataStatModule, self).end()

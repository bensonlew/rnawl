# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
# last_modify: 20180320
import os
import re
import shutil
from biocluster.module import Module
from biocluster.core.exceptions import OptionError
class CleanStatModule(Module):
    """
    原始数据统计，统计冗余度，fastq_dup.py(输入为文件夹时文件夹内需要有list.txt)
    统计碱基质量信息及Q20、Q30，fastx_quality_stats、q20q30_stat.py(输入为文件夹时文件夹内需要有list.txt)
    统计fastq序列基本信息，FastqStat.jar(输入为文件夹时文件夹内需要有list.txt)
    """

    def __init__(self, work_id):
        super(CleanStatModule, self).__init__(work_id)
        options = [
            {"name": "raw_list", "type": "infile", "format": "bacgenome.list_file"},  #rawdata list
            {"name": "clean_list", "type": "infile", "format": "bacgenome.list_file"},  #cleandata list
            {"name": "sample_name", "type": "string"},  # 样品名称
        ]
        self.add_option(options)
        self.step.add_steps("fastq_stat", "fastx")
        self.samples = {}
        self.moduels = []

    def set_step(self, event):
        if "start" in event["data"].keys():
            event["data"]["start"].start()
        if "end" in event["data"].keys():
            event["data"]["end"].finish()
        self.step.update()

    def check_options(self):
        if not self.option("raw_list").is_set:
            raise OptionError("必须设置进行原始的文件list", code="21401201")
        if not self.option("clean_list").is_set:
            raise OptionError("必须设置进行质控的文件list", code="21401202")
        if not self.option("sample_name"):
            raise OptionError("必须样品名称", code="21401203")
        self.option("clean_list").get_info()
        self.samples = self.option("clean_list").prop["samples"]

    def run_fastq_stat(self):
        options = {
            "raw_list": self.option("raw_list"),
            "clean_list": self.option("clean_list"),
            "sample_name": self.option("sample_name"),
        }
        self.fastq_stat = self.add_tool("bacgenome.bac_read_statistic")
        self.fastq_stat.set_options(options)
        self.moduels.append(self.fastq_stat)
        self.fastq_stat.on("start", self.set_step, {"start": self.step.fastq_stat})
        self.fastq_stat.on("end", self.set_step, {"end": self.step.fastq_stat})
        self.fastq_stat.on('end', self.set_output, 'fastq_stat')
        self.fastq_stat.run()

    def run_fastx(self):
        base_dir = os.path.dirname(self.option("clean_list").prop['path'])
        for s in self.samples.keys():
            if re.search('PE',s):
                options = {
                    "fastq": base_dir + '/' + self.samples[s][1]
                }
                self.fastx = self.add_tool("bacgenome.fastx_v2")
                self.fastx.set_options(options)
                self.moduels.append(self.fastx)
                self.fastx.on("start", self.set_step, {"start": self.step.fastx})
                self.fastx.on("end", self.set_step, {"end": self.step.fastx})
                self.fastx.on('end', self.set_output, 'fastx_{}_r'.format(s))
                self.fastx.run()
                options = {
                    "fastq": base_dir + '/' + self.samples[s][0]
                }
                self.fastx = self.add_tool("bacgenome.fastx_v2")
                self.fastx.set_options(options)
                self.moduels.append(self.fastx)
                self.fastx.on("start", self.set_step, {"start": self.step.fastx})
                self.fastx.on("end", self.set_step, {"end": self.step.fastx})
                self.fastx.on('end', self.set_output, 'fastx_{}_l'.format(s))
                self.fastx.run()

    def set_output(self, event):
        self.logger.info("开始set output")
        obj = event["bind_object"]
        m2 = re.match(r"fastx_(.+)", event["data"])
        if m2:
            new_dir = os.path.join(self.output_dir, "fastx")
            if not os.path.exists(new_dir):
                os.mkdir(new_dir)
            for f in os.listdir(obj.output_dir):
                n1 = re.match(r".+fastxstat$", f)
                if n1:
                    if os.path.exists(new_dir + "/" + m2.group(1) + ".clean_fastxstat"):
                        os.remove(new_dir + "/" + m2.group(1) + ".clean_fastxstat")
                    os.link(obj.output_dir + "/" + f, new_dir + "/" + m2.group(1) + ".clean_fastxstat")
        m3 = re.match(r"fastq_stat", event["data"])
        if m3:
            file =os.path.join(self.output_dir, 'data_QC')
            if not os.path.exists(file):
                os.mkdir(file)
            for f in os.listdir(obj.output_dir):
                n1 = re.match(r"qc_stat", f)
                if n1:
                    if os.path.exists(file + '/' + self.option('sample_name') + "_Illumina_statistics.xls"):
                        os.remove(file + '/' + self.option('sample_name') + "_Illumina_statistics.xls")
                    os.link(obj.output_dir + "/" + f, file + '/' + self.option('sample_name') + "_Illumina_statistics.xls")

    def run(self):
        super(CleanStatModule, self).run()
        self.run_fastq_stat()
        self.run_fastx()
        self.on_rely(self.moduels, self.end)

    def end(self):
        super(CleanStatModule, self).end()
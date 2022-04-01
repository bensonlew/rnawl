# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
# last_modify: 20210315
import os
import re
import shutil
from biocluster.module import Module
from biocluster.core.exceptions import OptionError
class CleanStatModule(Module):
    """
    真菌二代和三代统计
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
        self.fastq_stat = self.add_tool("bacgenome.bac_read_statistic")
        self.clean_fastx = self.add_module("fungi_genome.clean_fastx")
        self.samples = {}
        self.modules = []

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
        self.fastq_stat.set_options(options)
        self.fastq_stat.on('end', self.set_output, 'fastq_stat')
        self.fastq_stat.run()

    def run_fastx(self):
        options = {
            "clean_list": self.option("clean_list"),
        }
        self.clean_fastx.set_options(options)
        self.clean_fastx.on('end', self.set_output, 'fastx')
        self.clean_fastx.run()


    def set_output(self, event):
        self.logger.info("开始set output")
        obj = event["bind_object"]
        m2 = re.match(r"fastx", event["data"])
        if m2:
            new_dir = os.path.join(self.output_dir, "fastx")
            if not os.path.exists(new_dir):
                os.mkdir(new_dir)
            for f in os.listdir(obj.output_dir+"/fastx"):
                n1 = re.match(r".+clean_fastxstat$", f)
                if n1:
                    if os.path.exists(new_dir + "/" + f):
                        os.remove(new_dir + "/" + f)
                    os.link(obj.output_dir + "/fastx/" + f, new_dir + "/" + f)
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
        self.on_rely([self.clean_fastx, self.fastq_stat], self.end)

    def end(self):
        super(CleanStatModule, self).end()
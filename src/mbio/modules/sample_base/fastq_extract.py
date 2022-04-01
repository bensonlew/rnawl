#!/usr/bin/env python
# -*- coding: utf-8 -*-
# __author__ = "shijin"

from __future__ import division
import os
import shutil
from biocluster.module import Module


class FastqExtractModule(Module):
    """
    version 1.0
    author: shijin
    last_modify: 2017.9.25
    """
    def __init__(self, work_id):
        super(FastqExtractModule, self).__init__(work_id)
        options = [
            {"name": "in_fastq", "type": "infile", "format": "sequence.fastq,sequence.fastq_dir"},
            {"name": "output_fq", "type": "outfile", "format": "sequence.fastq_dir"},
            {"name": "output_length", "type": "outfile", "format": "sample.data_dir"},
            {"name": "output_list", "type": "outfile", "format": "sequence.info_txt"}

        ]
        self.add_option(options)
        # self.length_stat = self.add_tool("sample_base.length_stat")
        self.samples = []
        self.tools = []

    def check_options(self):
        return True

    def fastq_run(self):
        opts = {
            "in_fastq": self.option("in_fastq")
        }
        run_tool = self.add_tool("sample_base.fastq_extract")
        run_tool.set_options(opts)
        self.tools.append(run_tool)
        run_tool.on("end", self.move_to_output)
        run_tool.run()

    def fastq_dir_run(self):
        fq_list = []
        if self.option("in_fastq").format == "sequence.fastq_dir":
            self.samples = self.option("in_fastq").prop["fastq_basename"]
        for f in self.samples:
            self.logger.info("开始运行样本{}".format(f))
            fq_path = self.option("in_fastq").prop["path"] + "/" + f
            fq_list.append(fq_path)
        self.logger.info("开始进行样本合并")
        self.logger.info(" ".join(fq_list))
        os.system("cat {} > {}".format(" ".join(fq_list), self.work_dir + "/cat_sample.fq"))
        # change by wzy 20170925, join中的内容改为list
        self.logger.info("样本合并完成")
        opts = {
            "in_fastq": self.work_dir + "/cat_sample.fq"  # change by wzy 20170925, 输入文件改为合并后的文件
        }
        run_tool = self.add_tool("sample_base.fastq_extract")
        run_tool.set_options(opts)
        self.tools.append(run_tool)
        run_tool.on("end", self.move_to_output)
        run_tool.run()

    def run(self):
        # self.length_stat.on("end", self.end)  # change by wzy 20170925, 注释掉
        if self.option("in_fastq").format == "sequence.fastq":
            self.fastq_run()
        else:
            self.logger.info("输入文件为文件夹，开始进行并行运行")
            self.fastq_dir_run()
        super(FastqExtractModule, self).run()

    def end(self):
        # super(FastqExtractModule, self).end()
        self.option("output_fq").set_path(self.output_dir + "/fastq")
        self.option("output_length").set_path(self.output_dir + "/length")
        self.option("output_list").set_path(self.work_dir + "/info.txt")
        super(FastqExtractModule, self).end()

    def move_to_output(self):
        self.logger.info("进入移动文件过程")
        shutil.rmtree(self.output_dir)
        os.mkdir(self.output_dir)
        os.mkdir(self.output_dir + "/fastq")
        fq_dir = self.output_dir + "/fastq"
        os.mkdir(self.output_dir + "/length")
        length_dir = self.output_dir + "/length"
        if os.path.exists(self.work_dir + "/info.txt"):
            os.remove(self.work_dir + "/info.txt")
        with open(self.work_dir + "/info.txt", "a") as a:
            a.write("#file\tsample\tworkdir\tseqs_num\tbase_num\tmean_length\tmin_length\tmax_length\n")
        for tool in self.tools:
            for file in os.listdir(tool.output_dir + "/fastq"):
                file_path = os.path.join(tool.output_dir + "/fastq", file)
                file_name = os.path.split(file_path)[1]
                if not os.path.exists(fq_dir + "/" + file_name):
                    os.link(file_path, fq_dir + "/" + file_name)
                else:
                    with open(fq_dir + "/" + file_name, "a") as a:
                        content = open(file_path, "r").read()
                        a.write(content)
            for file in os.listdir(tool.output_dir + "/length"):
                file_path = os.path.join(tool.output_dir + "/length", file)
                file_name = os.path.split(file_path)[1]
                if not os.path.exists(length_dir + "/" + file_name):
                    os.link(file_path, length_dir + "/" + file_name)
                else:
                    with open(length_dir + "/" + file_name, "a") as a:
                        content = open(file_path, "r").read()
                        a.write(content)
            list_path = tool.work_dir + "/info.txt"  # 由于上面已默认将相同样本合并，此处应有调整
            with open(self.work_dir + "/info.txt", "a") as a:
                f = open(list_path, "r")
                f.readline()
                a.write(f.read())
        self.option('output_fq').set_path(self.output_dir + '/fastq')
        self.option('output_list').set_path(self.work_dir + '/info.txt')
        self.end()
        # self.get_length_stat()  # change by wzy 20170925, 以下注释掉

    # def get_length_stat(self):
    #     self.logger.info("移动样本文件过程结束，进入长度统计步骤")
    #     opts = {
    #         "length_dir": self.output_dir + "/length",
    #         "file_sample_list": self.work_dir + "/info.txt"  # 仅用来获取序列的最大长度
    #     }
    #     self.length_stat.set_options(opts)
    #     self.length_stat.on("end", self.end)
    #     self.length_stat.run()


#!/usr/bin/env python
# -*- coding: utf-8 -*-
# __author__ = "shijin"

from __future__ import division
import os,shutil,re
from biocluster.module import Module
from biocluster.core.exceptions import OptionError
#from biocluster.config import Config


class SampleExtractModule(Module):
    """
    version 1.0
    author: shijin
    last_modify: 2017.04.11
    """
    def __init__(self, work_id):
        super(SampleExtractModule, self).__init__(work_id)
        options = [
            {"name": "in_fastq", "type": "infile", "format": "sequence.fastq,sequence.fastq_dir"},
            {"name": "file_list", "type": "string", "default": "null"},
            {"name": "table_id", "type": "string", "default": ""},
            {"name": "info_txt", "type":"infile", "format": "sequence.info_txt"},
            {"name": "file_sample_list", "type": "outfile", "format": "sequence.info_txt"}
        ]
        self.add_option(options)
        self.samples = []
        self.tools = []
        self.info_path = []
        
    def check_options(self):
        return True
        
    def fastq_run(self, fastqfile):
        opts = {
            "in_fastq": fastqfile
        }
        run_tool = self.add_tool("meta.fastq_sample_extract")
        run_tool.set_options(opts)
        self.tools.append(run_tool)
        run_tool.on("end", self.end)
        run_tool.run()
        
    def fastq_dir_run(self):
        # if self.option("table_id") != "":
        self.samples = self.option("in_fastq").prop["fastq_basename"]
        for f in self.samples:
            fq_path = self.option("in_fastq").prop["path"] + "/" + f
            opts = {
                "in_fastq": fq_path
            }
            run_tool = self.add_tool("meta.fastq_sample_extract")
            run_tool.set_options(opts)
            self.tools.append(run_tool)
        if len(self.tools) >= 1:
            self.on_rely(self.tools, self.paste)
        elif len(self.tools) == 1:
            self.tools[0].on("end", self.paste)
        for tool in self.tools:
            tool.run()
        
    def run(self):
        super(SampleExtractModule, self).run()
        if self.option("in_fastq").format == "sequence.fastq":
            self.fastq_run(self.option("in_fastq"))
        else:
            self.logger.info("run fastq dir")
            self.fastq_dir_run()

    def paste(self):
        """
        for dir in os.listdir(self.work_dir):
            if dir.startswith("FastqSampleExtract"):
                dir_path = os.path.join(self.work_dir, dir + "/" + "info.txt")
                if os.path.isfile(dir_path):  # 增加判断，防止因为重运行导致没有这个文件 by ghd @ 20181024
                    self.info_path.append(dir_path)
                else:
                    self.logger.info("路径:%s,不存在，确定此任务是否有重运行" % dir_path)
        """
        if self.option("in_fastq").format == "sequence.fastq":
            dir_path = os.path.join(self.tools.work_dir, "info.txt")
            self.info_path.append(dir_path)
        else:
            for tool in self.tools:
                dir_path = os.path.join(tool.work_dir, "info.txt")
                self.info_path.append(dir_path)
        list_path = os.path.join(self.work_dir, "info.txt")
        self.logger.info("already here")
        with open(list_path, "w") as w:
            w.write("#file_path\tsample\twork_dir_path\tseq_num\tbase_num\tmean_length\tmin_length\tmax_length\n")
            for path in self.info_path:
                with open(path, "r") as r:
                    r.readline()
                    info_part = r.read()
                    w.write(info_part)
        self.end()

    def end(self):
        if self.option("file_list") == "null":
            if not os.path.exists(self.work_dir + "/info.txt"):
                # 此时输入文件为sequence.fastq self.tools只有一个tool by ghd @ 20181021
                os.link(os.path.join(self.tools[0].work_dir, "info.txt"), self.work_dir + "/info.txt")
                # os.link(self.work_dir + "/FastqSampleExtract/info.txt", self.work_dir + "/info.txt")
        #if self.option("file_list") == "null" and self.option("table_id") != "":
        #    self.logger.info(self.option("table_id"))
        #    self.set_sample_db()
        #    self.option("file_sample_list").set_path(Config().WORK_DIR + "/sample_data/" +
        #                                             self.option("table_id") + "/info.txt")
        #else:
        self.option("file_sample_list").set_path(self.work_dir + "/info.txt")
        with open(self.option("file_sample_list").prop["path"], "r") as file:
            file.readline()
            try:
                next(file)
            except:
                raise OptionError("样本检测没有找到样本，请重新检查文件的改名信息", code="")
        super(SampleExtractModule, self).end()

    '''
    def set_sample_db(self):
        os.mkdir(Config().WORK_DIR + "/sample_data/" + self.option("table_id"))
        table_dir = os.path.join(Config().WORK_DIR + "/sample_data", self.option("table_id"))
        new_info_path = os.path.join(table_dir, "info.txt")
        old_info_path = self.work_dir + "/info.txt"
        with open(new_info_path, "w") as w:
            with open(old_info_path, "r") as r:
                w.write("#file_path\tsample\twork_dir_path\tseq_num\tbase_num\tmean_length\tmin_length\tmax_length\n")
                r.readline()
                for line in r:
                    line = line.strip()
                    lst = line.split("\t")
                    sample_name = lst[1]
                    """
                    file_name = os.path.basename(lst[0])
                    sample_name = lst[1]
                    key = file_name + "::" + sample_name
                    if key in file_list.keys():
                    """
                    new_tool_lst = lst[2].split("/")
                    new_tool_path = table_dir + "/" + new_tool_lst[-1]
                    self.mv(lst[2], new_tool_path, sample_name)
                    w.write(lst[0] + "\t" + sample_name + "\t" + new_tool_path + "\t" + lst[3] + "\t" + lst[
                        4] + "\t" + lst[5] + "\t" + lst[6] + "\t" + lst[7] + "\n")
    '''

    """
    def mv(self, old_path, new_path, key):
        if self.option("file_list") != "null":
            file_list = eval(self.option("file_list"))
            old_name = file_list[key][1]
            new_name = file_list[key][0]
            if new_name.find(".") != -1 or new_name.find(" ") != -1:
                raise Exception("样本名称中带.和空格，请更改样本名称为不带.的名称后再进行流程")
            if new_name == "OUT" or new_name == "IN":
                raise Exception("样本名称不能为IN与OUT")
        else:
            old_name = key
            new_name = key
        if not os.path.exists(new_path):
            os.mkdir(new_path)
        output_path = new_path + "/output"
        if not os.path.exists(output_path):
            os.mkdir(output_path)
        output_length_path = output_path + "/length"
        output_fa_path = output_path + "/fa"
        if not os.path.exists(output_length_path):
            os.mkdir(output_length_path)
        if not os.path.exists(output_fa_path):
            os.mkdir(output_fa_path)
        for file in os.listdir(old_path + "/output/length"):
            sample = re.match("(.+)\.length_file", file).group(1)
            if sample == old_name:
                old_file = os.path.join(old_path + "/output/length", file)
                new_file = os.path.join(output_length_path, new_name + ".length_file")
                if os.path.exists(new_file):
                    with open(new_file, "a") as a:
                        with open(old_file, "r") as r:
                            for line in r:
                                a.write(line)
                else:
                    shutil.copy(old_file, new_file)  # 此处不能用os.link()
        for file in os.listdir(old_path + "/output/fa"):
            sample = re.match("(.+)\.fa", file).group(1)
            if sample == old_name:
                old_file = os.path.join(old_path + "/output/fa", file)
                new_file = os.path.join(output_fa_path, new_name + ".fasta")
                if os.path.exists(new_file):
                    with open(new_file, "a") as a:
                        with open(old_file, "r") as r:
                            for line in r:
                                a.write(line)
                else:
                    shutil.copy(old_file, new_file)
    """
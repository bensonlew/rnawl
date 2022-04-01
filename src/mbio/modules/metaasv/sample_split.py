#!/usr/bin/env python
# -*- coding: utf-8 -*-
# __author__ = "qingchen.zhang"

from __future__ import division
import os
from biocluster.module import Module
from mbio.packages.metaasv.common_function import link_dir, link_file


class SampleSplitModule(Module):
    """
    metaasv 对fastq序列进行拆分
    """
    def __init__(self, work_id):
        super(SampleSplitModule, self).__init__(work_id)
        options = [
            {"name": "in_fastq", "type": "infile", "format": "sequence.fastq,sequence.fastq_dir"},
            # {"name": "info_txt", "type":"infile", "format": "sequence.info_txt"},
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
        self.run_tool = self.add_tool("metaasv.split_fastq")
        self.run_tool.set_options(opts)
        self.run_tool.on("end", self.set_output)
        self.run_tool.run()
        
    def fastq_dir_run(self):
        self.samples = self.option("in_fastq").prop["fastq_basename"]
        for f in self.samples:
            fq_path = self.option("in_fastq").prop["path"] + "/" + f
            opts = {
                "in_fastq": fq_path
            }
            run_tool = self.add_tool("metaasv.split_fastq")
            run_tool.set_options(opts)
            self.tools.append(run_tool)
        if len(self.tools) >= 1:
            self.on_rely(self.tools, self.set_output)
        elif len(self.tools) == 1:
            self.tools[0].on("end", self.set_output)
        for tool in self.tools:
            tool.run()
        
    def run(self):
        super(SampleSplitModule, self).run()
        if self.option("in_fastq").format == "sequence.fastq":
            self.fastq_run(self.option("in_fastq"))
        else:
            self.logger.info("run fastq dir")
            self.fastq_dir_run()

    def set_output(self):
        if len(os.listdir(self.output_dir)) == 0:
            if self.option("in_fastq").format == "sequence.fastq_dir":
                ww = open(self.work_dir +"/list.txt", "w")
                n = 1
                for tool in self.tools:
                    for file in os.listdir(tool.output_dir):
                        file_path = os.path.join(tool.output_dir, file)
                        ## 为防止多批次的样本合并出现问题，下面增加代码记录改名前和改名后的名称
                        ## list第一列为改名前，第二列为改名后
                        new_file_path = os.path.join(self.output_dir, file)
                        origin_path = tool.option("in_fastq").prop['path']
                        if os.path.exists(new_file_path):
                            new_file_path = os.path.join(self.output_dir, "duplicate"+ str(n) +file)
                            ww.write("duplicate"+ str(n)+file+ "\t" +file+ "\t"+ origin_path+"\n")
                            n += 1
                        else:
                            ww.write(file + "\t" + file +"\t"+ origin_path+ "\n")
                        link_file(file_path, new_file_path)

            else:
                ww = open(self.work_dir +"/list.txt", "w")
                n = 1
                for file in os.listdir(self.run_tool.output_dir):
                    file_path = os.path.join(self.run_tool.output_dir, file)
                    new_file_path = os.path.join(self.output_dir, file)
                    origin_path = self.run_tool.option("in_fastq").prop['path']
                    if os.path.exists(new_file_path):
                        new_file_path = os.path.join(self.output_dir, "duplicate"+ str(n) +file)
                        ww.write("duplicate"+ str(n)+file+ "\t" +file+ "\t"+ origin_path+"\n")
                        n += 1
                    else:
                        ww.write(file + "\t" + file +"\t"+ origin_path+ "\n")
                    link_file(file_path, new_file_path)
        self.end()

    def end(self):
        super(SampleSplitModule, self).end()
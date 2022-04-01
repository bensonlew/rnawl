# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# last_modify: 20171210

import os
from biocluster.module import Module
from biocluster.core.exceptions import OptionError


class SingleNcrnaQcModule(Module):
    """
    进行单个ncRNA的fastq的质控：是否去除前三个碱基，保留前51个，用cutadapt进行质控
    """
    def __init__(self, work_id):
        super(SingleNcrnaQcModule, self).__init__(work_id)
        options = [
            {"name": "fastq", "type": "infile", "format": "sequence.fastq"},  # fastq文件
            {"name": "fastq_l", "type": "infile", "format": "sequence.fastq"},  # 样本左端fastq文件
            {"name": "fastq_r", "type": "infile", "format": "sequence.fastq"},  # 样本右端fastq文件
            {"name": "min_length", "type": "string", "default": "18"},  # reads最短长度设置，舍弃长度小于此值的序列
            {"name": "low_quality_base", "type": "string", "default": "33"},  # 过滤掉序列两端的低质量碱基数
            {"name": "r_adaptor", "type": "string", "default": "TGGAATTCTCGGGTGCCAAGG"},  # 左端接头
            {"name": "l_adaptor", "type": "string", "default": "GATCGTCGGACTGTAGAACTCTGAAC"},  # 右端接头
            {"name": "adaptor", "type": "string", "default": "GATCGTCGGACTGTAGAACTCTGAAC"},  # 单端接头
            {"name": "cut_left", "type": "string", "default": "False"}  # miRNA是否要切除前3bp，保留前51bp
        ]
        self.add_option(options)

    def finish_update(self, event):
        step = getattr(self.step, event['data'])
        step.finish()
        self.step.update()

    def check_options(self):
        if self.option("fastq").is_set:
            pass
        elif self.option("fastq_l").is_set and self.option("fastq_r").is_set:
            pass
        else:
            raise OptionError("请设置fastq文件或者fastq_l和fastq_r文件")

    def run_cut_left(self):
        if self.option("fastq").is_set:
            options = {
                "fastq": self.option("fastq")
            }
        else:
            options = {
                "fastq_l": self.option("fastq_l"),
                "fastq_r": self.option("fastq_r")
            }
        self.cut_left = self.add_tool("datasplit.cut_left51")
        self.cut_left.set_options(options)
        self.step.add_steps("cut_left")
        step = getattr(self.step, "cut_left")
        step.start()
        self.step.update()
        self.cut_left.on('end', self.finish_update, "cut_left")
        self.cut_left.on("end", self.run_cutadapt)
        self.cut_left.run()

    def run_cutadapt(self):
        options = {
            "min_length": self.option("min_length"),
            "low_quality_base": self.option("low_quality_base"),
            "r_adaptor": self.option("r_adaptor"),
            "l_adaptor": self.option("l_adaptor"),
            "adaptor": self.option("adaptor")
        }
        if self.option("fastq").is_set:
            if self.option("cut_left") is "True":
                options["fastq"] = self.cut_left.option("cut_fastq")
            else:
                options["fastq"] = self.option("fastq")
        else:
            if self.option("cut_left") is "True":
                options["fastq_l"] = self.cut_left.option("cut_fastq_l")
                options["fastq_r"] = self.cut_left.option("cut_fastq_r")
            else:
                options["fastq_l"] = self.option("fastq_l")
                options["fastq_r"] = self.option("fastq_r")
        self.cutadapt = self.add_tool("datasplit.cutadapt")
        self.cutadapt.set_options(options)
        self.step.add_steps("cutadapt")
        step = getattr(self.step, "cutadapt")
        step.start()
        self.step.update()
        self.cutadapt.on('end', self.finish_update, "cutadapt")
        self.cutadapt.on('end', self.set_output, 'cutadapt')
        self.cutadapt.run()

    def set_output(self, event):
        obj = event["bind_object"]
        if event["data"] == "cutadapt":
            if not os.path.isdir(self.cutadapt.output_dir):
                raise Exception("需要移动到output目录的文件夹不存在")
            for f in os.listdir(self.cutadapt.output_dir):
                os.link(os.path.join(self.cutadapt.output_dir, f), os.path.join(self.output_dir, f))
            self.end()

    def run(self):
        super(SingleNcrnaQcModule, self).run()
        if self.option("cut_left") is "True":
            self.run_cut_left()
        else:
            self.run_cutadapt()

    def end(self):
        super(SingleNcrnaQcModule, self).end()

# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# last_modify: 20171210

import os
import re
from biocluster.module import Module
from biocluster.core.exceptions import OptionError


class NcrnaQcModule(Module):
    """
    进行多个ncRNA的fastq的质控：并行投递module：single_ncrna_qc
    """
    def __init__(self, work_id):
        super(NcrnaQcModule, self).__init__(work_id)
        options = [
            {"name": "list_file", "type": "infile", "format": "datasplit.list_file"},  # fastq文件及对应的样本
            {"name": "fq_type", "type": "string", "default": "PE"},  # 质控方式，PE/SE
            {"name": "minlen", "type": "string", "default": "18"},  # reads最短长度设置，舍弃长度小于此值的序列
            {"name": "low_quality_base", "type": "string", "default": "33"},  # 过滤掉序列两端的低质量碱基数
            {"name": "adaptor", "type": "string", "default": "GATCGTCGGACTGTAGAACTCTGAAC"},  # 单端接头
            {"name": "r_adaptor", "type": "string", "default": "TGGAATTCTCGGGTGCCAAGG"},  # 右端接头
            {"name": "l_adaptor", "type": "string", "default": "GATCGTCGGACTGTAGAACTCTGAAC"},  # 左端接头
            {"name": "cut_left", "type": "string", "default": "False"}  # miRNA是否要切除前3bp，保留前51bp
        ]
        self.add_option(options)
        self.step.add_steps("single_ncrna_qc")
        self.fastqs = {}
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
            raise OptionError("必须设置fastq文件的list")
        self.fastqs = self.option("list_file").prop["samples"]
        if self.option("fq_type") not in ["SE", "PE"]:
            raise OptionError("fq_type只能是SE/PE")
        if self.option("fq_type") == "SE":
            for s in self.fastqs.keys():
                if len(self.fastqs[s]) == 2:
                    raise OptionError("list_file中的类型是PE，和fq_type不一致，请检查")
        else:
            for s in self.fastqs.keys():
                if len(self.fastqs[s]) == 1:
                    raise OptionError("list_file中的类型是SE，和fq_type不一致，请检查")

    def run_ncrna_qc(self):
        options = {
            "min_length": self.option("minlen"),
            "low_quality_base": self.option("low_quality_base"),
            "adaptor": self.option("adaptor"),
            "r_adaptor": self.option("r_adaptor"),
            "l_adaptor": self.option("l_adaptor"),
            "cut_left": self.option("cut_left")
        }
        modules = []
        for s in self.fastqs.keys():
            if self.option("fq_type") == "SE":
                options["fastq"] = self.fastqs[s][0]
            else:
                options["fastq_l"] = self.fastqs[s][0]
                options["fastq_r"] = self.fastqs[s][1]
            self.single_ncrna_qc = self.add_module("datasplit.single_ncrna_qc")
            self.single_ncrna_qc.set_options(options)
            self.single_ncrna_qc.on('start', self.set_step, {'start': self.step.single_ncrna_qc})
            self.single_ncrna_qc.on('end', self.set_step, {'end': self.step.single_ncrna_qc})
            self.single_ncrna_qc.on('end', self.set_output, 'ncrna_qc_{}'.format(s))
            modules.append(self.single_ncrna_qc)
            self.start_times += 1
        for m in modules:
            m.run()

    def set_output(self, event):
        obj = event["bind_object"]
        self.end_times += 1
        sample = event["data"].split("ncrna_qc_")[1]
        if not os.path.isdir(obj.output_dir):
            raise Exception("需要移动到output目录的文件夹不存在")
        for f in os.listdir(obj.output_dir):
            m = re.search(r"sample(_.+).fastq", f)
            if m:
                f1 = sample + m.group(1) + ".fastq"
            else:
                raise Exception("输出clean fastq命名不符合要求")
            os.link(os.path.join(obj.output_dir, f), os.path.join(self.output_dir, f1))
        if self.start_times == self.end_times:
            self.end()

    def run(self):
        super(NcrnaQcModule, self).run()
        self.run_ncrna_qc()

    def end(self):
        super(NcrnaQcModule, self).end()

# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# modified 2018.04.25

from biocluster.core.exceptions import OptionError
from biocluster.module import Module
import os


class SsrPrimerModule(Module):
    """
    样本基因组SSR分析及引物设计
    """
    def __init__(self, work_id):
        super(SsrPrimerModule, self).__init__(work_id)
        options = [
            {"name": "fastq_list", "type": "string"},  # 样本对应的fastq
            {"name": "tm1", "type": "float", "default": 57.0},  # float, Tm1 (℃)
            {"name": "tm2", "type": "float", "default": 63.0},  # float, Tm2 (℃),要大于tm1
            {"name": "product_size", "type": "string", "default": "300-500"},  # Product Size(bp),数值要为int,范围间用-分隔，多个用,分隔,如：100-300
            {"name": "primer_num", "type": "int", "default": 3},  # Max Pairs Primer Number,范围:[1,5]
        ]
        self.add_option(options)
        self.start_times, self.end_times = 0, 0

    def check_options(self):
        if not self.option("fastq_list"):
            raise OptionError("请设置fastq_list", code="24501301")

    def run_make_config(self):
        """
        tool:make_config
        """
        options = {
            "fastq_list": self.option("fastq_list")
        }
        self.make_config = self.add_tool("wgs.make_config")
        self.make_config.set_options(options)
        self.make_config.on("end", self.run_soap_denovo)
        self.make_config.run()

    def run_soap_denovo(self):
        """
        tool:soap_denovo
        """
        self.logger.info(self.make_config.output_dir)
        for f in os.listdir(self.make_config.output_dir):
            self.start_times += 1
            config_file = os.path.join(self.make_config.output_dir, f)
            sample_name = f.split(".config")[0]
            options = {
                "config_file": config_file
            }
            self.soap_denovo = self.add_tool("wgs.soap_denovo")
            self.soap_denovo.set_options(options)
            self.soap_denovo.on("end", self.run_ssr_misa, sample_name)
            self.soap_denovo.on("end", self.set_output, "soap_" + sample_name)
            self.soap_denovo.run()

    def run_ssr_misa(self, event):
        """
        tool:ssr_misa
        """
        obj = event["bind_object"]
        sample_name = event["data"]
        scafseq_file = os.path.join(obj.output_dir, sample_name + ".denovo.scafSeq")
        options = {
            "scafseq_file": scafseq_file
        }
        self.ssr_misa = self.add_tool("wgs.ssr_misa")
        self.ssr_misa.set_options(options)
        self.ssr_misa.on("end", self.run_ssr_stat, sample_name)
        self.ssr_misa.on("end", self.run_ssr_primer_design, sample_name)
        self.ssr_misa.run()

    def run_ssr_stat(self, event):
        """
        tool:ssr_stat
        """
        obj = event["bind_object"]
        sample_name = event["data"]
        misa_file = os.path.join(obj.output_dir, sample_name + ".denovo.scafSeq.misa")
        options = {
            "misa_file": misa_file
        }
        self.ssr_stat = self.add_tool("wgs.ssr_stat")
        self.ssr_stat.set_options(options)
        self.ssr_stat.on("end", self.set_output, "ssr_stat")
        self.ssr_stat.run()

    def run_ssr_primer_design(self, event):
        """
        tool:ssr_primer_design
        """
        obj = event["bind_object"]
        sample_name = event["data"]
        misa_file = os.path.join(obj.output_dir, sample_name + ".denovo.scafSeq.misa")
        scafseq_file = os.path.join(obj.output_dir, sample_name + ".denovo.scafSeq")
        options = {
            "misa_file": misa_file,
            "scafseq_file": scafseq_file,
            "tm1": self.option("tm1"),
            "tm2": self.option("tm2"),
            "product_size": self.option("product_size"),
            "primer_num": self.option("primer_num")
        }
        self.ssr_primer_design = self.add_tool("wgs.ssr_primer_design")
        self.ssr_primer_design.set_options(options)
        self.ssr_primer_design.on("end", self.set_output, "ssr_primer")
        self.ssr_primer_design.run()

    def set_output(self, event):
        obj = event["bind_object"]
        self.end_times += 1
        ssr_compare = os.path.join(self.output_dir, "ssr_compare")
        if not os.path.exists(ssr_compare):
            os.mkdir(ssr_compare)
        if event["data"].startswith("soap_"):
            sample = event["data"].split("soap_")[1]
            os.link(os.path.join(obj.output_dir, sample + ".denovo.scafSeq"), os.path.join(ssr_compare, sample + ".denovo.scafSeq"))
        else:
            for f in os.listdir(obj.output_dir):
                if f.endswith(".result"):
                    os.link(os.path.join(obj.output_dir, f), os.path.join(ssr_compare, f))
                else:
                    os.link(os.path.join(obj.output_dir, f), os.path.join(self.output_dir, f))
        if self.start_times * 3 == self.end_times:
            self.end()

    def run(self):
        self.run_make_config()
        super(SsrPrimerModule, self).run()

    def end(self):
        super(SsrPrimerModule, self).end()

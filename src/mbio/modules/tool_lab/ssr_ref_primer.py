# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# modified 2018.04.25

from biocluster.core.exceptions import OptionError
from biocluster.module import Module
import os


class SsrRefPrimerModule(Module):
    """
    样本基因组SSR分析及引物设计
    """
    def __init__(self, work_id):
        super(SsrRefPrimerModule, self).__init__(work_id)
        options = [
            {"name": "reffa", "type": "infile", "format": "sequence.fasta"},  # 参考基因组ref.fa
            {"name": "tm1", "type": "float", "default": 57.0},  # float, Tm1 (℃)
            {"name": "tm2", "type": "float", "default": 63.0},  # float, Tm2 (℃),要大于tm1
            {"name": "product_size", "type": "string", "default": "300-500"},  # Product Size(bp),数值要为int,范围间用-分隔，多个用,分隔,如：100-300
            {"name": "primer_num", "type": "int", "default": 3},  # Max Pairs Primer Number,范围:[1,5]
            {"name": "needini", "type": "bool", "default": False},  # 做misa的时候是否要配置ini文件，默认不配
            {"name": "needprimer", "type": "bool", "default": True},  # 是否要做引物，默认做
            {"name": "rept_1", "type": "int", "default": 10},
            {"name": "rept_2", "type": "int", "default": 6},
            {"name": "rept_3", "type": "int", "default": 5},
            {"name": "rept_4", "type": "int", "default": 5},
            {"name": "rept_5", "type": "int", "default": 5},
            {"name": "rept_6", "type": "int", "default": 5},
            {"name": "ssr_distance", "type": "int", "default": 100},
        ]
        self.add_option(options)
        self.start_times, self.end_times = 0, 0
        self.ssr_ref_primer_design, self.ssr_misa, self.ssr_stat = '', '', ''

    def check_options(self):
        if not self.option("reffa").is_set:
            raise OptionError("请设置reffa", code="24501401")

    # def run_make_config(self):
    #     """
    #     tool:make_config
    #     """
    #     options = {
    #         "reffa": self.option("reffa")
    #     }
    #     self.make_config = self.add_tool("wgs.make_config")
    #     self.make_config.set_options(options)
    #     self.make_config.on("end", self.run_soap_denovo)
    #     self.make_config.run()

    # def run_soap_denovo(self):
    #     """
    #     tool:soap_denovo
    #     """
    #     self.logger.info(self.make_config.output_dir)
    #     for f in os.listdir(self.make_config.output_dir):
    #         self.start_times += 1
    #         config_file = os.path.join(self.make_config.output_dir, f)
    #         sample_name = f.split(".config")[0]
    #         options = {
    #             "config_file": config_file
    #         }
    #         self.soap_denovo = self.add_tool("wgs.soap_denovo")
    #         self.soap_denovo.set_options(options)
    #         self.soap_denovo.on("end", self.run_ssr_misa, sample_name)
    #         self.soap_denovo.run()

    def run_ssr_misa(self):
        """
        tool:ssr_misa 生成misa文件；
        """
        options = {
            "scafseq_file": self.option("reffa"),
            "needini": self.option("needini"),
            "rept_1": self.option("rept_1"),
            "rept_2": self.option("rept_2"),
            "rept_3": self.option("rept_3"),
            "rept_4": self.option("rept_4"),
            "rept_5": self.option("rept_5"),
            "rept_6": self.option("rept_6"),
            "ssr_distance": self.option("ssr_distance"),
        }
        self.ssr_misa = self.add_tool("tool_lab.ssr_misa")
        self.ssr_misa.set_options(options)
        self.ssr_misa.on("end", self.run_ssr_stat)
        if self.option("needprimer"):
            self.ssr_misa.on("end", self.run_ssr_ref_primer_design)  # ssr_ref_primer_design.pl
        self.ssr_misa.on("end", self.set_output, "ssr_misa")
        self.ssr_misa.run()

    def run_ssr_stat(self, event):
        """
        tool:ssr_stat

        """
        obj = event["bind_object"]
        sample_name = event["data"]
        # misa_file = os.path.join(obj.output_dir, sample_name + ".misa")
        misa_file = os.path.join(obj.output_dir, os.path.basename(self.option("reffa").prop["path"]) + ".misa")
        options = {
            "misa_file": misa_file
        }
        self.ssr_stat = self.add_tool("wgs.ssr_stat")
        self.ssr_stat.set_options(options)
        self.ssr_stat.on("end", self.set_output, "ssr_stat")
        self.ssr_stat.run()

    def run_ssr_ref_primer_design(self, event):
        """
        tool:ssr_ref_primer_design
        """
        obj = event["bind_object"]
        sample_name = event["data"]
        # misa_file = os.path.join(obj.output_dir, sample_name + ".denovo.scafSeq.misa")
        # scafseq_file = os.path.join(obj.output_dir, sample_name + ".denovo.scafSeq")
        misa_file = os.path.join(obj.output_dir, os.path.basename(self.option("reffa").prop["path"]) + ".misa")
        options = {
            "misa_file": misa_file,
            "scafseq_file": self.option("reffa"),
            "tm1": self.option("tm1"),
            "tm2": self.option("tm2"),
            "product_size": self.option("product_size"),
            "primer_num": self.option("primer_num")
        }
        self.ssr_ref_primer_design = self.add_tool("wgs.ssr_ref_primer_design")
        self.ssr_ref_primer_design.set_options(options)
        self.ssr_ref_primer_design.on("end", self.set_output, "ssr_primer")
        self.ssr_ref_primer_design.run()

    def set_output(self, event):
        obj = event["bind_object"]
        self.end_times += 1
        for f in os.listdir(obj.output_dir):
            if f == 'ref.fa':
                continue
            if os.path.exists(os.path.join(self.output_dir, f)):
                os.remove(os.path.join(self.output_dir, f))
            os.link(os.path.join(obj.output_dir, f), os.path.join(self.output_dir, f))
        if self.option("needprimer"):
            if self.end_times == 3:
                self.end()
        else:
            if self.end_times == 2:
                self.end()

    def run(self):
        self.run_ssr_misa()
        super(SsrRefPrimerModule, self).run()

    def end(self):
        super(SsrRefPrimerModule, self).end()

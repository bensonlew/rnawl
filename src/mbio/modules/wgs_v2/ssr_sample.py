# -*- coding: utf-8 -*-
# __author__: zengjing
# modified: 20190321

import os
from biocluster.module import Module
from biocluster.core.exceptions import OptionError


class SsrSampleModule(Module):
    """
    对单个样本进行SSR设计及统计
    """
    def __init__(self, work_id):
        super(SsrSampleModule, self).__init__(work_id)
        options = [
            {"name": "bam_file", "type": "infile", "format": "align.bwa.sam", "required": True},  # 样本对应的bam文件
            # {"name": "ref_chrlist", "type": "infile", "format": "wgs_v2.bcf", "required": True},  # cankao基因组对应的ref.chrlist
            {"name": "ref_chrlist", "type": "string", "required": True},  # cankao基因组对应的ref.chrlist
            # {"name": "ref_misa", "type": "infile", "format": "wgs_v2.bcf", "required": True},  # 参考基因组的SSR:ref.fa.misa
            {"name": "ref_misa", "type": "string", "required": True},  # 参考基因组的SSR:ref.fa.misa
            {"name": "sample_name", "type": "string"}
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("bam_file").is_set:
            raise OptionError("请设置样本的bam文件")
        # if not self.option("ref_chrlist").is_set:
        #     raise OptionError("请设置参考基因组的ref.chrlist")
        # if not self.option("ref_misa").is_set:
        #     raise OptionError("请设置参考基因组的ref.fa.misa")

    def run_bam_split(self):
        options = {
            "bam_file": self.option("bam_file"),
            "ref_chrlist": self.option("ref_chrlist")
        }
        self.bam_split = self.add_tool("wgs_v2.bam_split")
        self.bam_split.set_options(options)
        self.bam_split.on("end", self.run_bam_ssr)
        self.bam_split.run()

    def run_bam_ssr(self):
        self.all_ssr = []
        for f in os.listdir(self.bam_split.output_dir):
            bam_file = os.path.join(self.bam_split.output_dir, f)
            options = {
                "bam_file": bam_file,
                "ref_misa": self.option("ref_misa")
            }
            self.bam_ssr = self.add_tool("wgs_v2.bam_ssr")
            self.bam_ssr.set_options(options)
            self.all_ssr.append(self.bam_ssr)
        if len(self.all_ssr) == 1:
            self.all_ssr[0].on("end", self.run_ssr_misa)
        else:
            self.on_rely(self.all_ssr, self.run_ssr_misa)
        for bam_ssr in self.all_ssr:
            bam_ssr.run()

    def run_ssr_misa(self):
        ssr_list = os.path.join(self.work_dir, "ssr.list")
        with open(ssr_list, "wb") as w:
            for bam_ssr in self.all_ssr:
                for f in os.listdir(bam_ssr.output_dir):
                    if f.endswith("ssr.result"):
                        w.write(f.split(".ssr.result")[0] + "\t" + os.path.join(bam_ssr.output_dir, f) + "\n")
        options = {
            "ssr_list": ssr_list,
            "sample_name": self.option("sample_name")
        }
        self.ssr_misa = self.add_tool("wgs_v2.ssr_misa")
        self.ssr_misa.set_options(options)
        self.ssr_misa.on("end", self.set_output, "ssr_misa")
        self.ssr_misa.run()

    def set_output(self, event):
        obj = event["bind_object"]
        for f in os.listdir(obj.output_dir):
            if os.path.exists(os.path.join(self.output_dir, f)):
                os.remove(os.path.join(self.output_dir, f))
            os.link(os.path.join(obj.output_dir, f), os.path.join(self.output_dir, f))
        self.end()

    def run(self):
        super(SsrSampleModule, self).run()
        sample_name = self.option("sample_name")
        self.run_bam_split()

    def end(self):
        self.add_upload_dir(self.output_dir)
        super(SsrSampleModule, self).end()

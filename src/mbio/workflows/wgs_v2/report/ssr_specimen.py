# -*- coding: utf-8 -*-
# __author__: zengjing
# modified: 20190321

import os
from bson import ObjectId
from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError


class SsrSpecimenWorkflow(Workflow):
    """
    交互分析：所有样本的SSR设计
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(SsrSpecimenWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "bam_list", "type": "infile", "format": "wgs_v2.bam_list"},  # 所有样本的bam.list
            {"name": "bam_path", "type": "infile", "format": "wgs_v2.bam_dir"},  # 所有样本的bam文件夹路径
            # {"name": "ref_chrlist", "type": "infile", "format": "wgs_v2.bcf", "required": True},  # cankao基因组对应的ref.chrlist
            {"name": "ref_chrlist", "type": "string", "required": True},  # cankao基因组对应的ref.chrlist
            # {"name": "ref_misa", "type": "infile", "format": "wgs_v2.bcf", "required": True},  # 参考基因组的SSR:ref.fa.misa
            {"name": "ref_misa", "type": "string", "required": True},  # 参考基因组的SSR:ref.fa.misa
            {"name": "update_info", "type": "string"},
            {"name": "main_id", "type": "string"},
            {"name": "task_id", "type": "string"}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())

    def check_options(self):
        if not self.option("bam_list").is_set and not self.option("bam_path").is_set:
            raise OptionError("请设置bam.list或者bam_path")
        # if not self.option("ref_chrlist").is_set:
        #     raise OptionError("请设置参考基因组的ref.chrlist")
        # if not self.option("ref_misa").is_set:
        #     raise OptionError("请设置参考基因组的ref.fa.misa")

    def get_bam_list(self):
        """
        得到bam.list
        """
        if self.option("bam_list").is_set:
            self.bam_list = self.option("bam_list").prop["path"]
        else:
            self.bam_list = os.path.join(self.work_dir, "bam.list")
            with open(self.bam_list, "wb") as w:
                for f in os.listdir(self.option("bam_path").prop["path"]):
                    if f.endswith("sort.bam"):
                        s = f.split(".sort.bam")[0]
                        w.write(s + "\t" + os.path.join(self.option("bam_path").prop["path"], f) + "\n")

    def run_ssr_sample(self):
        self.all_ssr = []
        with open(self.bam_list, "rb") as f:
            for line in f:
                item = line.strip().split("\t")
                options = {
                    "bam_file": item[1],
                    "ref_chrlist": self.option("ref_chrlist"),
                    "ref_misa": self.option("ref_misa"),
                    "sample_name": item[0]
                }
                self.ssr_sample = self.add_module("wgs_v2.ssr_sample")
                self.ssr_sample.set_options(options)
                self.all_ssr.append(self.ssr_sample)
        if len(self.all_ssr) == 1:
            self.all_ssr[0].on("end", self.all_ssr_stat)
        else:
            self.on_rely(self.all_ssr, self.all_ssr_stat)
        for sample_ssr in self.all_ssr:
            sample_ssr.run()

    def all_ssr_stat(self):
        self.ssr_list = os.path.join(self.work_dir, "ssr_result.list")
        all_ssr = os.path.join(self.output_dir, "all.ssr.stat.xls")
        with open(self.ssr_list, "wb") as w, open(all_ssr, "wb") as w1:
            w1.write("Sample ID\tSSR Number\tc\tc*\tp1\tp2\tp3\tp4\tp5\tp6\n")
            for sample_ssr in self.all_ssr:
                for f in os.listdir(sample_ssr.output_dir):
                    if f.endswith("ssr.result"):
                        w.write(f.split(".ssr.result")[0] + "\t" + os.path.join(sample_ssr.output_dir, f) + "\n")
                    elif f.endswith("ssr.stat.xls"):
                        with open(os.path.join(sample_ssr.output_dir, f), "rb") as f1:
                            lines = f1.readlines()
                            w1.write(lines[1])
        self.run_ssr_merge()

    def run_ssr_merge(self):
        options = {
            "ref_misa": self.option("ref_misa"),
            "ssr_list": self.ssr_list
        }
        self.ssr_merge = self.add_tool("wgs_v2.ssr_merge")
        self.ssr_merge.set_options(options)
        self.ssr_merge.on("end", self.set_output, "ssr_merge")
        self.ssr_merge.run()

    def set_output(self, event):
        for f in os.listdir(event["bind_object"].output_dir):
            new = os.path.join(self.output_dir, f)
            if os.path.exists(new):
                os.remove(new)
            os.link(os.path.join(event["bind_object"].output_dir, f), new)
        self.set_db()

    def set_db(self):
        self.logger.info("开始进行导表")
        if self.option("main_id"):
            ssr_api = self.api.api("wgs_v2.ssr_specimen")
            ssr_stat_path = os.path.join(self.output_dir, "all.ssr.stat.xls")
            ssr_vcf = os.path.join(self._sheet.output, "final.ssr.result.vcf")
            ssr_api.add_sg_ssr_marker_stat_specimen(task_id=self.option("task_id"), ssr_id=self.option("main_id"), ssr_stat_path=ssr_stat_path)
            ssr_api.update_info(coll="sg_ssr_marker", query_dict={"main_id": ObjectId(self.option("main_id"))}, update_dict={"ssr_vcf": ssr_vcf})
        self.end()

    def run(self):
        self.get_bam_list()
        self.run_ssr_sample()
        super(SsrSpecimenWorkflow, self).run()

    def end(self):
        self.add_upload_dir(self.output_dir)
        super(SsrSpecimenWorkflow, self).end()

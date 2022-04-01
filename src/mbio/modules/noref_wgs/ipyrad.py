# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# last_modify: 20181225

import os
from biocluster.module import Module
from biocluster.core.exceptions import OptionError


class IpyradModule(Module):
    """
    ipyrad：
    """
    def __init__(self, work_id):
        super(IpyradModule, self).__init__(work_id)
        options = [
            {"name": "fastq_dir", "type": "infile", "format": "bsa.dir", "required": True},
            {"name": "enzyme_method", "type": "string", "required": True},  # 酶切方案,GBS/RAD
            {"name": "cutseq", "type": "string", "required": True},  # 酶切序列, 如：TCGA, TTAA
        ]
        self.add_option(options)
        self.end_times = 0

    def check_options(self):
        if self.option("enzyme_method") not in ["RAD", "GBS"]:
            raise OptionError("酶切方案: %s只能是RAD/GBS" , variables=( self.option("enzyme_method")), code="25500205")
        if self.option("enzyme_method") == "GBS":
            if len(self.opiton("cutseq").split(", ")) != 2:
                raise OptionError("酶切方案为GBS的时候cutseq:%s 必须为2个" , variables=( self.option("enzyme_method")), code="25500206")

    def run_ipyrad(self):
        options = {
            "fastq_dir": self.option("fastq_dir"),
            "analysis_step": "1234567",
            "enzyme_method": self.option("enzyme_method"),
            "cutseq": self.option("cutseq")
        }
        self.ipyrad = self.add_tool("noref_wgs.ipyrad")
        self.ipyrad.set_options(options)
        self.ipyrad.on("end", self.run_tag_stat)
        self.ipyrad.on("end", self.run_snp_stat)
        self.ipyrad.on("end", self.set_output, "ipyrad")
        self.ipyrad.run()

    def run_tag_stat(self):
        cluster_dir = self.ipyrad.output_dir + "/data_clust_0.85"
        data_loci = self.ipyrad.output_dir + "/data_outfiles/data.loci"
        options = {
            "cluster_dir": cluster_dir,
            "data_loci": data_loci,
            "total_sample_num": self.total_sample_num
        }
        self.tag_stat = self.add_tool("noref_wgs.ipyrad_tag_stat")
        self.tag_stat.set_options(options)
        self.tag_stat.on("end", self.set_output, "tag_stat")
        self.tag_stat.run()

    def run_snp_stat(self):
        vcf_path = self.ipyrad.output_dir + "/data_outfiles/data.vcf"
        options = {
            "vcf_path": vcf_path
        }
        self.snp_stat = self.add_tool("noref_wgs.snp_stat")
        self.snp_stat.set_options(options)
        self.snp_stat.on("end", self.set_output, "snp_stat")
        self.snp_stat.run()

    def link_file(self, old, new):
        if os.path.exists(new):
            os.remove(new)
        os.link(old, new)

    def set_output(self, event):
        self.logger.info("设置结果目录")
        obj = event["bind_object"]
        if event["data"] == "ipyrad":
            data_vcf = os.path.join(obj.output_dir, "data_outfiles/data.vcf")
            self.link_file(data_vcf, os.path.join(self.output_dir, "data.vcf"))
        else:
            for f in os.listdir(obj.output_dir):
                self.link_file(os.path.join(obj.output_dir, f), os.path.join(self.output_dir, f))
        self.end_times += 1
        if self.end_times == 3:
            self.end()

    def run(self):
        super(IpyradModule, self).run()
        fq_list = os.path.join(self.option("fastq_dir").prop["path"], "fq.list")
        if os.path.exists(fq_list):
            self.total_sample_num = len(open(fq_list, 'r').readlines())
        elif self.option("enzyme_method") == "GBS":
            self.total_sample_num = len(glob.glob(r".*fastq.gz")) / 2
        else:
            self.total_sample_num = len(glob.glob(r".*fastq.gz"))
        self.run_ipyrad()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        super(IpyradModule, self).end()

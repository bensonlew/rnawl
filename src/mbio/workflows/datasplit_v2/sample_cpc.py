# -*- coding: utf-8 -*-
# __author__ = "zengjing"
# last_modify: 20191012

import os
import json
import gevent
from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError


class SampleCpcWorkflow(Workflow):
    """
    对拆分板里的rna(常规RNA)/dna(动植物基因组)/prokaryotic_rna(原核mRNA)/lncrna(LncRNA)/microbial_genome(微生物基因组质控)样本进行CPC
    CPC：随机抽取10000条序列，比对到nt看样本的物种污染情况，比对到rfam，看样本的核糖体污染情况
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(SampleCpcWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "sample_list", "type": "infile", "format": "datasplit.list_file"},  # 样本名称、样本路径，双端用分号分隔，分号前是l，后是r端
            {"name": "split_id", "type": "string"},
            {"name": "update_info", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.end_times = 0

    def check_options(self):
        if not self.option("sample_list").is_set:
            raise OptionError("请设置sample_list文件")

    def get_all_specimen(self):
        """
        得到每个样本对应的序列文件
        """
        self.specimen_path = {}
        self.nt_specimen, self.rfam_specimen = [], []
        with open(self.option("sample_list").prop["path"], "rb") as f:
            for line in f:
                item = line.strip().split("\t")
                self.specimen_path[item[0]] = item[1].split(";")
                if item[0] in self.nt_specimen:
                    raise Exception("样本:%s重复，请检查" % item[0])
                self.nt_specimen.append(item[0])
                if item[2] == "microbial_genome":
                    continue
                if item[0] in self.rfam_specimen:
                    raise Exception("样本:%s重复，请检查" % item[0])
                self.rfam_specimen.append(item[0])

    def run_blast_nt(self):
        """
        将样本比对到nt数据库，进行统计，看样本的污染情况
        """
        for sample in self.nt_specimen:
            if len(self.specimen_path[sample]) == 1:
                options = {"fastq": self.specimen_path[sample][0]}
            elif(len(self.specimen_path[sample]) == 2):
                options = {
                    "fastq_l": self.specimen_path[sample][0],
                    "fastq_r": self.specimen_path[sample][1],
                }
            # options["database"] = "nt"
            # options["database"] = "nt_v20200604"
            options["database"] = "nt_v20211009"
            # options["database"] = "nt_v202107"
            self.blast_nt = self.add_module("datasplit_v2.fastq_blast")
            self.blast_nt.set_options(options)
            self.blast_nt.on("end", self.set_output, sample)
            self.blast_nt.run()

    def run_blast_rfam(self):
        """
        将样本比对到rfam数据库，看样本的核糖体污染情况
        """
        for sample in self.rfam_specimen:
            if len(self.specimen_path[sample]) == 1:
                options = {
                    "fastq": self.specimen_path[sample][0],
                    "evalue": 1e-5,
                    "num": "100000",
                    "num_alignment": 1,
                    }
            elif(len(self.specimen_path[sample]) == 2):
                options = {
                    "fastq_l": self.specimen_path[sample][0],
                    "fastq_r": self.specimen_path[sample][1],
                    "evalue": 1e-5,
                    "num": "100000",
                    "num_alignment": 1,
                }
            # options["database"] = "rfam"
            options["database"] = "rfam_v14.6"
            self.blast_rfam = self.add_module("datasplit_v2.fastq_blast")
            self.blast_rfam.set_options(options)
            self.blast_rfam.on("end", self.set_output, sample)
            self.blast_rfam.run()

    def set_output(self, event):
        """
        输出结果目录整理
        """
        obj = event["bind_object"]
        sample_name = event["data"]
        for f in os.listdir(obj.output_dir):
            f1 = os.path.join(obj.output_dir, f)
            f2 = os.path.join(self.output_dir, sample_name + "--" + f)
            if os.path.exists(f2):
                os.remove(f2)
            os.link(f1, f2)
        self.end_times += 1
        if self.end_times == len(self.nt_specimen)+len(self.rfam_specimen):
            self.set_db()

    def set_db(self):
        """
        保存结果到mongo数据库
        """
        self.logger.info("保存结果到mongo数据库")
        if self.option("split_id"):
            datasplit_api = self.api.api("datasplit.datasplit_new")
            datasplit_api.update_cpc_info(self.option("split_id"), self.output_dir)
        self.end()

    def run(self):
        self.get_all_specimen()
        self.run_blast_nt()
        self.run_blast_rfam()
        if len(self.rfam_specimen) + len(self.nt_specimen) == 0:
            self.start_listener()
            self.fire("start")
            gevent.spawn_later(5, self.end)
            self.end()
        super(SampleCpcWorkflow, self).run()

    def end(self):
        # self.add_upload_dir(self.output_dir)
        super(SampleCpcWorkflow, self).end()
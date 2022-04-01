# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# modified 2018.08.29

from biocluster.core.exceptions import OptionError
from biocluster.module import Module
import re
import os


class SingleMappingModule(Module):
    """
    单个样本的mapping比对,注：ref.fa要好建索引，单个样本可能有多个批次的，在这儿会将所有的样本合并为一个bam
    """
    def __init__(self, work_id):
        super(SingleMappingModule, self).__init__(work_id)
        options = [
            {"name": "fastq_list", "type": "infile", "format": "wgs.fastq_list"},  #
            # fastq路径list.txt文件，第一列样本名称，第二列左端序列，第三列右端序列，第四列文库类型
            {"name": "ref_fa", "type": "infile", "format": "sequence.fasta"},
            {"name": "specimen_id", "type": "string"},  # 样本名称
            {"name": "num", "type": "string"},   # 每个样本的bam or sam文件中的rg id 必须唯一
            {"name": "large_genome", "type": "string", "default": "false"}
        ]
        self.add_option(options)
        self.samtools_view_tools, self.bwa_mem_tools = [], []

    def check_options(self):
        if not self.option("fastq_list").is_set:
            raise OptionError("请设置fastq_list", code="24502019")
        if not self.option("specimen_id"):
            raise OptionError("请设置样本名称", code="24502020")
        if not self.option("ref_fa").is_set:
            raise OptionError("请设置参考序列", code="24502021")
        else:
            ref_dir = os.path.dirname(self.option("ref_fa").prop["path"])
            base_name = os.path.basename(self.option("ref_fa").prop["path"])
            amb_path = os.path.join(ref_dir, base_name + ".amb")
            ann_path = os.path.join(ref_dir, base_name + ".ann")
            bwt_path = os.path.join(ref_dir, base_name + ".bwt")
            fai_path = os.path.join(ref_dir, base_name + ".fai")
            pac_path = os.path.join(ref_dir, base_name + ".pac")
            sa_path = os.path.join(ref_dir, base_name + ".sa")
            if not os.path.exists(amb_path):
                raise OptionError("ref_fa的索引不全，缺少ref.fa.amb", code="24502022")
            if not os.path.exists(ann_path):
                raise OptionError("ref_fa的索引不全，缺少ref.fa.ann", code="24502023")
            if not os.path.exists(bwt_path):
                raise OptionError("ref_fa的索引不全，缺少ref.fa.bwt", code="24502024")
            if not os.path.exists(fai_path):
                raise OptionError("ref_fa的索引不全，缺少ref.fa.fai", code="24502025")
            if not os.path.exists(pac_path):
                raise OptionError("ref_fa的索引不全，缺少ref.fa.pac", code="24502026")
            if not os.path.exists(sa_path):
                raise OptionError("ref_fa的索引不全，缺少ref.fa.sa", code="24502027")

    def get_sample_list(self):
        """
        筛选出WGS/RAD文库，非RAD文库要跑tool：picard_mkdup
        """
        self.specimen_ids = {}
        self.specimen_list = []
        self.specimen_type = ""
        with open(self.option("fastq_list").prop["path"], "r") as f:
            lines = f.readlines()
            for line in lines:
                item = line.strip().split("\t")
                self.specimen_type = item[3]
                self.specimen_list.append(item[0])
                self.specimen_ids[item[0]] = {"l": item[1], "r": item[2]}
        self.logger.info(self.specimen_ids)

    def run_bwa_mem(self):
        """
        tool:bwa_mem
        """
        for s in self.specimen_list:
            options = {
                "fastq_l": self.specimen_ids[s]["l"],
                "fastq_r": self.specimen_ids[s]["r"],
                "sample_name": s,
                "ref_fa": self.option("ref_fa"),
                "num": self.option('num')
            }
            self.bwa_mem = self.add_tool("wgs.bwa_mem")
            self.bwa_mem.set_options(options)
            self.bwa_mem.on("end", self.run_samtools_view, s)
            self.bwa_mem.run()

    def run_samtools_view(self, event):
        """
        tool: samtools_view
        """
        obj = event["bind_object"]
        for f in os.listdir(obj.output_dir):
            sam_file = os.path.join(obj.output_dir, f)
        options = {
            "sam_file": sam_file
        }
        self.samtools_view = self.add_tool("wgs.samtools_view")
        self.samtools_view.set_options(options)
        self.samtools_view.on("end", self.set_output, "view")
        self.samtools_view_tools.append(self.samtools_view)
        if len(self.samtools_view_tools) == len(self.specimen_list):
            if len(self.samtools_view_tools) == 1:
                self.samtools_view_tools[0].on("end", self.run_samtools_merge, "single")
            else:
                self.on_rely(self.samtools_view_tools, self.run_samtools_merge, "more")
            for tool in self.samtools_view_tools:
                tool.run()

    def run_samtools_merge(self, event):
        """
        tool: samtools_merge
        """
        view_list = []
        if event["data"] == "single":
            for f in os.listdir(self.samtools_view.output_dir):
                view_list.append(os.path.join(self.samtools_view.output_dir, f))
        else:
            view_dir = os.path.join(self.work_dir, "samtools_view")
            for f in os.listdir(view_dir):
                view_list.append(os.path.join(view_dir, f))
        options = {
            "bam_list": ';'.join(view_list),
            "specimen_id": self.option("specimen_id")
        }
        self.samtools_merge = self.add_tool("wgs.samtools_merge")
        self.samtools_merge.set_options(options)
        self.samtools_merge.on("end", self.run_samtools_sort)
        self.samtools_merge.run()

    def run_samtools_sort(self):
        """
        tool: samtools_sort
        """
        for f in os.listdir(self.samtools_merge.output_dir):
            merged_bam_file = os.path.join(self.samtools_merge.output_dir, f)
        options = {
            "merged_bam_file": merged_bam_file
        }
        self.samtools_sort = self.add_tool("wgs.samtools_sort")
        self.samtools_sort.set_options(options)
        self.samtools_sort.on("end", self.run_samtools_index)
        self.samtools_sort.run()

    def run_samtools_index(self):
        """
        tool: samtools_index
        """
        for f in os.listdir(self.samtools_sort.output_dir):
            sort_bam_file = os.path.join(self.samtools_sort.output_dir, f)
        options = {
            "sort_bam_file": sort_bam_file,
            "ref_fa": self.option("ref_fa"),
        }
        self.samtools_index = self.add_tool("wgs.samtools_index")
        self.samtools_index.set_options(options)
        if self.specimen_type == "WGS":
            self.samtools_index.on("end", self.run_picard_mkdup)
        else:
            self.samtools_index.on("end", self.set_output, "out_bam")
        self.samtools_index.on("end", self.set_output, "sort_bams")
        self.samtools_index.run()

    def run_picard_mkdup(self):
        """
        tool: picard_mkdup
        """
        for f in os.listdir(self.samtools_index.output_dir):
            if f.endswith(".bam"):
                sort_bam_file = os.path.join(self.samtools_index.output_dir, f)
        options = {
            "sort_bam_file": sort_bam_file,
            "mkdup_method": "picard" if self.option("large_genome") == 'false' else 'samtools'
        }
        self.picard_mkdup = self.add_tool("wgs.picard_mkdup")
        self.picard_mkdup.set_options(options)
        self.picard_mkdup.on("end", self.set_output, "out_bam")
        self.picard_mkdup.run()

    def linkdir(self, old, new):
        if os.path.isfile(old):
            if os.path.exists(new):
                os.remove(new)
            os.link(old, new)
        else:
            if not os.path.exists(new):
                os.mkdir(new)
            for f in os.listdir(old):
                f1 = os.path.join(old, f)
                f1_ = os.path.join(new, f)
                self.linkdir(f1, f1_)

    def set_output(self, event):
        """
        """
        obj = event["bind_object"]
        if event["data"] == "view":
            view_dir = os.path.join(self.work_dir, "samtools_view")
            if not os.path.exists(view_dir):
                os.mkdir(view_dir)
            self.linkdir(obj.output_dir, view_dir)
        elif event["data"] == "sort_bams":
            sort_path = os.path.join(self.output_dir, "sort_bams")
            if not os.path.exists(sort_path):
                os.mkdir(sort_path)
            self.linkdir(obj.output_dir, sort_path)
        else:
            self.linkdir(obj.output_dir, self.output_dir)
        if event["data"] == "out_bam":
            self.end()

    def run(self):
        super(SingleMappingModule, self).run()
        self.get_sample_list()
        self.run_bwa_mem()

    def end(self):
        super(SingleMappingModule, self).end()

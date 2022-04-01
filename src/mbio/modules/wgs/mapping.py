# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# modified 2018.04.08

from biocluster.core.exceptions import OptionError
from biocluster.module import Module
import re
import os


class MappingModule(Module):
    """
    mapping比对,注：ref.fa要好建索引
    """
    def __init__(self, work_id):
        super(MappingModule, self).__init__(work_id)
        options = [
            {"name": "fastq_list", "type": "infile", "format": "wgs.fastq_list"},
            # fastq路径list.txt文件，第一列样本名称，第二列左端序列，第三列右端序列，第四列文库类型
            {"name": "ref_fa", "type": "infile", "format": "sequence.fasta"},
            {"name": "large_genome", "type": "string", "default": "false"}
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("fastq_list").is_set:
            raise OptionError("请设置fastq_list", code="24500901")
        if not self.option("ref_fa").is_set:
            raise OptionError("请设置参考序列", code="24500902")
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
                raise OptionError("ref_fa的索引不全，缺少ref.fa.amb", code="24500903")
            if not os.path.exists(ann_path):
                raise OptionError("ref_fa的索引不全，缺少ref.fa.ann", code="24500904")
            if not os.path.exists(bwt_path):
                raise OptionError("ref_fa的索引不全，缺少ref.fa.bwt", code="24500905")
            if not os.path.exists(fai_path):
                raise OptionError("ref_fa的索引不全，缺少ref.fa.fai", code="24500906")
            if not os.path.exists(pac_path):
                raise OptionError("ref_fa的索引不全，缺少ref.fa.pac", code="24500907")
            if not os.path.exists(sa_path):
                raise OptionError("ref_fa的索引不全，缺少ref.fa.sa", code="24500908")

    def get_sample_list(self):
        """
        筛选出WGS/RAD文库，非RAD文库要跑tool：picard_mkdup
        """
        self.specimen_ids, self.specimen_path = {}, {}
        self.specimen_list = []
        with open(self.option("fastq_list").prop["path"], "r") as f:
            lines = f.readlines()
            for line in lines:
                item = line.strip().split("\t")
                specimen_id = item[0].split("-")[0]
                if specimen_id not in self.specimen_ids:
                    self.specimen_ids[specimen_id] = {}
                    self.specimen_list.append(specimen_id)
                self.specimen_ids[specimen_id][item[0]] = {"l": item[1], "r": item[2], "type": item[3]}
        self.logger.info(self.specimen_ids)
        self.logger.info(self.specimen_list)
        specimen_dir = os.path.join(self.work_dir, "specimen_dir")
        if not os.path.exists(specimen_dir):
            os.mkdir(specimen_dir)
        for s in self.specimen_list:
            s_path = os.path.join(specimen_dir, s + ".list")
            with open(s_path, "w") as w:
                for s1 in self.specimen_ids[s].keys():
                    w.write(s1 + "\t" + self.specimen_ids[s][s1]["l"] + "\t" + self.specimen_ids[s][s1]["r"] + "\t" + self.specimen_ids[s][s1]["type"] + "\n")
            self.specimen_path[s] = s_path

    def run_single_mapping(self):
        """
        single_mapping:一个分析样本进行mapping
        """
        num = 1
        self.single_mapping_modules = []
        specimen_module = open(self.work_dir + "/specimen_module_id.list", "a")
        for s in self.specimen_list:
            options = {
                "fastq_list": self.specimen_path[s],
                "ref_fa": self.option("ref_fa"),
                "specimen_id": s,
                "num": str(num),
                "large_genome": self.option('large_genome')
            }
            num += 1
            self.single_mapping = self.add_module("wgs.single_mapping")
            self.single_mapping.set_options(options)
            self.single_mapping.on("end", self.set_output, s)
            self.single_mapping_modules.append(self.single_mapping)
            specimen_module.write(s + "\t" + self.single_mapping.id + "\n")
        specimen_module.close()
        self.on_rely(self.single_mapping_modules, self.get_bam_list)
        for module in self.single_mapping_modules:
            module.run()

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
        self.linkdir(obj.output_dir, self.output_dir)

    def get_bam_list(self):
        with open(self.output_dir + "/bam.list", "w") as w:
            for f in os.listdir(self.output_dir):
                if f.endswith(".bam"):
                    sample_name = f.split(".")[0]
                    w.write(sample_name + "\t" + os.path.join(self.output_dir, f) + "\n")
        self.end()

    def run(self):
        super(MappingModule, self).run()
        self.logger.info("rrr:{}".format(self.option('large_genome')))
        if os.path.exists(self.output_dir + "/bam.list"):
            os.remove(self.output_dir + "/bam.list")
        self.get_sample_list()
        self.run_single_mapping()

    def end(self):
        super(MappingModule, self).end()

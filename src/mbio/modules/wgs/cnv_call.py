# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# modified 2018.04.11

from biocluster.core.exceptions import OptionError
from biocluster.module import Module
import os
import time


class CnvCallModule(Module):
    """
    call cnv
    """
    def __init__(self, work_id):
        super(CnvCallModule, self).__init__(work_id)
        options = [
            {"name": "bam_list", "type": "string"},  # bam文件路径list.txt文件，第一列样本名，第二列bam路径
            {"name": "ref_gff", "type": "string"},  # ref.gff
        ]
        self.add_option(options)
        self.start_times, self.end_times = 0, 0

    def check_options(self):
        if not self.option("bam_list"):
            raise OptionError("请设置bam list", code="24500501")
        if not self.option("ref_gff"):
            raise OptionError("请设置ref.gff文件", code="24500502")

    def get_bam_list(self):
        """
        得到bam_list
        """
        self.bam_list = {}
        self.specimen = {}
        self.all_specimen = []
        with open(self.option("bam_list"), "r") as f:
            for line in f:
                item = line.strip().split("\t")
                self.all_specimen.append(item[0])
                self.bam_list[item[0]] = item[1]
                name = item[0].split("-")[0]
                if name not in self.specimen.keys():
                    self.specimen[name] = []
                self.specimen[name].append(item[0])
        self.start_times = len(self.bam_list.keys()) * 3

    def run_cnvnator(self):
        """
        tool:cnvnator
        """
        for s in self.all_specimen:
            options = {
                "bam_file": self.bam_list[s]
            }
            self.cnvnator = self.add_tool("wgs.cnvnator")
            self.cnvnator.set_options(options)
            self.cnvnator.on("end", self.run_cnv_anno, s)
            self.cnvnator.on("end", self.set_output, s)
            self.cnvnator.run()

    def run_cnv_anno(self, event):
        """
        tool: cnv_anno
        """
        obj = event["bind_object"]
        for f in os.listdir(obj.output_dir):
            cnv_file = os.path.join(obj.output_dir, f)
        options = {
            "cnv_file": cnv_file,
            "ref_gff": self.option("ref_gff")
        }
        self.cnv_anno = self.add_tool("wgs.cnv_anno")
        self.cnv_anno.set_options(options)
        self.cnv_anno.on("end", self.run_cnv_stat, event["data"])
        self.cnv_anno.on("end", self.set_output, event["data"])
        self.cnv_anno.run()

    def run_cnv_stat(self, event):
        """
        tool: cnv_stat
        """
        obj = event["bind_object"]
        for f in os.listdir(obj.output_dir):
            cnv_anno = os.path.join(obj.output_dir, f)
        options = {
            "cnv_anno": cnv_anno
        }
        self.cnv_stat = self.add_tool("wgs.cnv_stat")
        self.cnv_stat.set_options(options)
        self.cnv_stat.on("end", self.set_output, event["data"])
        self.cnv_stat.run()

    def set_output(self, event):
        """
        """
        obj = event["bind_object"]
        length_path = os.path.join(self.output_dir, "length")
        anno_path = os.path.join(self.output_dir, "anno")
        cnv_path = os.path.join(self.output_dir, "cnv")
        if not os.path.exists(length_path):
            os.mkdir(length_path)
        if not os.path.exists(anno_path):
            os.mkdir(anno_path)
        if not os.path.exists(cnv_path):
            os.mkdir(cnv_path)
        for f in os.listdir(obj.output_dir):
            if f.endswith("length.xls"):
                if os.path.exists(os.path.join(length_path, f)):
                    os.remove(os.path.join(length_path, f))
                os.link(os.path.join(obj.output_dir, f), os.path.join(length_path, f))
            elif f.endswith("anno.xls"):
                if os.path.exists(os.path.join(anno_path, f)):
                    os.remove(os.path.join(anno_path, f))
                os.link(os.path.join(obj.output_dir, f), os.path.join(anno_path, f))
            elif f.endswith("cnv.xls"):
                if os.path.exists(os.path.join(cnv_path, f)):
                    os.remove(os.path.join(cnv_path, f))
                os.link(os.path.join(obj.output_dir, f), os.path.join(cnv_path, f))
            else:
                if os.path.exists(os.path.join(self.output_dir, f)):
                    os.remove(os.path.join(self.output_dir, f))
                os.link(os.path.join(obj.output_dir, f), os.path.join(self.output_dir, f))
        time.sleep(1)
        self.end_times += 1
        if self.start_times == self.end_times:
            self.get_total_cnv_stat()
            self.end()

    def get_total_cnv_stat(self):
        """
        得到所有样本的统计表
        """
        total_stat = os.path.join(self.output_dir, "cnv.stat.xls")
        with open(total_stat, "w") as w:
            header = "#Sample ID\tdeletion\tduplication\tgene\n"
            w.write(header)
            for s in self.specimen.keys():
                dele_num, dup_num, gene_num = 0, 0, 0
                for s1 in self.specimen[s]:
                    f = os.path.join(self.output_dir, s1 + ".cnv.stat.xls")
                    with open(f, "r") as r:
                        lines = r.readlines()
                        for line in lines[1:]:
                            item = line.strip().split("\t")
                            if item[0].startswith("del"):
                                dele_num += int(item[1])
                                gene_num += int(item[2])
                            if item[0].startswith("dup"):
                                dup_num += int(item[1])
                                gene_num += int(item[2])
                    os.remove(f)
                w.write(s1 + "\t" + str(dele_num) + "\t" + str(dup_num) + "\t" + str(gene_num) + "\n")

    def run(self):
        super(CnvCallModule, self).run()
        self.get_bam_list()
        self.run_cnvnator()

    def end(self):
        super(CnvCallModule, self).end()

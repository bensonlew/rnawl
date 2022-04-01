# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# modified 2018.04.13

from biocluster.core.exceptions import OptionError
from biocluster.module import Module
import os
import time


class SvCallModule(Module):
    """
    call sv
    """
    def __init__(self, work_id):
        super(SvCallModule, self).__init__(work_id)
        options = [
            {"name": "bam_list", "type": "string"},  # bam文件路径list.txt文件，第一列样本名，第二列bam路径
            {"name": "ref_gff", "type": "string"},
        ]
        self.add_option(options)
        self.start_times, self.end_times = 0, 0

    def check_options(self):
        if not self.option("bam_list"):
            raise OptionError("请设置bam list", code="24501501")
        if not self.option("ref_gff"):
            raise OptionError("请设置ref_gff", code="24501502")

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

    def run_breakdancer(self):
        """
        tool: breakdancer
        """
        for s in self.all_specimen:
            options = {
                "bam_file": self.bam_list[s]
            }
            self.breakdancer = self.add_tool("wgs.breakdancer")
            self.breakdancer.set_options(options)
            self.breakdancer.on("end", self.run_sv_anno, s)
            self.breakdancer.on("end", self.set_output, s)
            self.breakdancer.run()

    def run_sv_anno(self, event):
        """
        tool: sv_anno
        """
        obj = event["bind_object"]
        for f in os.listdir(obj.output_dir):
            sv_file = os.path.join(obj.output_dir, f)
        options = {
            "sv_file": sv_file,
            "ref_gff": self.option("ref_gff")
        }
        self.sv_anno = self.add_tool("wgs.sv_anno")
        self.sv_anno.set_options(options)
        self.sv_anno.on("end", self.run_sv_stat, event["data"])
        self.sv_anno.on("end", self.set_output, event["data"])
        self.sv_anno.run()

    def run_sv_stat(self, event):
        """
        tool: sv_stat
        """
        obj = event["bind_object"]
        for f in os.listdir(obj.output_dir):
            sv_anno = os.path.join(obj.output_dir, f)
        options = {
            "sv_anno": sv_anno
        }
        self.sv_stat = self.add_tool("wgs.sv_stat")
        self.sv_stat.set_options(options)
        self.sv_stat.on("end", self.set_output, event["data"])
        self.sv_stat.run()

    def set_output(self, event):
        """
        """
        obj = event["bind_object"]
        length_path = os.path.join(self.output_dir, "length")
        anno_path = os.path.join(self.output_dir, "anno")
        sv_path = os.path.join(self.output_dir, "sv")
        if not os.path.exists(length_path):
            os.mkdir(length_path)
        if not os.path.exists(anno_path):
            os.mkdir(anno_path)
        if not os.path.exists(sv_path):
            os.mkdir(sv_path)
        for f in os.listdir(obj.output_dir):
            if f.endswith("length.xls"):
                if os.path.exists(os.path.join(length_path, f)):
                    os.remove(os.path.join(length_path, f))
                os.link(os.path.join(obj.output_dir, f), os.path.join(length_path, f))
            elif f.endswith("anno.xls"):
                if os.path.exists(os.path.join(anno_path, f)):
                    os.remove(os.path.join(anno_path, f))
                os.link(os.path.join(obj.output_dir, f), os.path.join(anno_path, f))
            elif f.endswith("sv.xls"):
                if os.path.exists(os.path.join(sv_path, f)):
                    os.remove(os.path.join(sv_path, f))
                os.link(os.path.join(obj.output_dir, f), os.path.join(sv_path, f))
            else:
                if os.path.exists(os.path.join(self.output_dir, f)):
                    os.remove(os.path.join(self.output_dir, f))
                os.link(os.path.join(obj.output_dir, f), os.path.join(self.output_dir, f))
        time.sleep(1)
        self.end_times += 1
        if self.start_times == self.end_times:
            self.get_total_sv_stat()
            self.end()

    def get_total_sv_stat(self):
        """
        得到所有样本的统计表
        """
        total_stat = os.path.join(self.output_dir, "sv.stat.xls")
        with open(total_stat, "w") as w:
            header = "#Sample ID\tCTX\tDEL\tINS\tINV\tITX\tgene\n"
            w.write(header)
            for s in self.specimen.keys():
                ctx_num, del_num, ins_num, inv_num, itx_num, gene_num = 0, 0, 0, 0, 0, 0
                for s1 in self.specimen[s]:
                    f = os.path.join(self.output_dir, s1 + ".sv.stat.xls")
                    with open(f, "r") as r:
                        lines = r.readlines()
                        for line in lines[1:]:
                            item = line.strip().split("\t")
                            if item[0] == "CTX":
                                ctx_num += int(item[1])
                                gene_num += int(item[2])
                            if item[0] == "DEL":
                                del_num += int(item[1])
                                gene_num += int(item[2])
                            if item[0] == "INS":
                                ins_num += int(item[1])
                                gene_num += int(item[2])
                            if item[0] == "INV":
                                inv_num += int(item[1])
                                gene_num += int(item[2])
                            if item[0] == "ITX":
                                itx_num += int(item[1])
                                gene_num += int(item[2])
                    os.remove(f)
                w.write(s1 + "\t" + str(ctx_num) + "\t" + str(del_num) + "\t" + str(ins_num) + "\t" + str(inv_num) +\
                        "\t" + str(itx_num) + "\t" + str(gene_num) + "\n")

    def run(self):
        super(SvCallModule, self).run()
        self.get_bam_list()
        self.run_breakdancer()

    def end(self):
        super(SvCallModule, self).end()

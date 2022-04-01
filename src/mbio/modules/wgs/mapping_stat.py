# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# modified 2018.04.09

from biocluster.core.exceptions import OptionError
from biocluster.module import Module
import re
import os
import time


class MappingStatModule(Module):
    """
    mapping比对的统计,生成插入片段、测序深度、覆盖度图的数据
    注：ref.fa要好建索引
    """
    def __init__(self, work_id):
        super(MappingStatModule, self).__init__(work_id)
        options = [
            {"name": "bam_list", "type": "string"},  # bam文件路径list.txt文件，第一列样本名，第二列bam路径,若有mkdup_metric文件，则放到相同bam路径下
            {"name": "ref_dict", "type": "string"},  # ref.ref_dict
            {"name": "step_num", "type": "int", "default": 10000},
            {"name": "ref_fa", "type": "infile", "format": "sequence.fasta"},
            {"name": "large_genome", "type": "string", "default": "false"}
        ]
        self.add_option(options)
        self.start_times, self.end_times = 0, 0

    def check_options(self):
        if not self.option("bam_list"):
            raise OptionError("请设置bam_list", code="24501001")
        if not self.option("ref_dict"):
            raise OptionError("请设置ref_dict", code="24501002")

    def get_bam_list(self):
        """
        得到bam_list
        """
        self.bam_list = {}
        self.all_specimen = []
        with open(self.option("bam_list"), "r") as f:
            for line in f:
                item = line.strip().split("\t")
                self.all_specimen.append(item[0])
                self.bam_list[item[0]] = {"bam": item[1]}
                mkdup_metric = os.path.join(os.path.dirname(item[1]), item[0] + ".metric")
                if os.path.exists(mkdup_metric):
                    self.bam_list[item[0]]["metric"] = mkdup_metric
        self.start_times = len(self.bam_list.keys()) * 2
        self.logger.info(self.all_specimen)

    def run_samtools_stats(self):
        """
        tool: samtools_stats，得到
        """
        for s in self.all_specimen:
            options = {
                "bam_file": self.bam_list[s]["bam"],
                "ref_fa": self.option("ref_fa"),
                "large_genome": self.option('large_genome')
            }
            self.samtools_stats = self.add_tool("wgs.samtools_stats")
            self.samtools_stats.set_options(options)
            self.samtools_stats.on("end", self.run_stats_stat, s)
            self.samtools_stats.run()

    def run_stats_stat(self, event):
        """
        tool: coverage_stat, 对samtools_stats和picard_mkdup的结果进行统计
        """
        obj = event["bind_object"]
        sample_name = event["data"]
        for f in os.listdir(obj.output_dir):
            map_stat = os.path.join(obj.output_dir, f)
        options = {
            "map_stat": map_stat,
            "ref_dict": self.option("ref_dict")
        }
        if "metric" in self.bam_list[sample_name].keys():
            options["mkdup_metric"] = self.bam_list[sample_name]["metric"]
        self.map_stat = self.add_tool("wgs.coverage_stat")
        self.map_stat.set_options(options)
        self.map_stat.on("end", self.set_output)
        self.map_stat.run()

    def run_samtools_depth(self):
        """
        tool: samtools_depth, 覆盖度数据
        """
        for s in self.bam_list.keys():
            options = {
                "bam_file": self.bam_list[s]["bam"],
                "ref_fa": self.option("ref_fa"),
                "large_genome": self.option('large_genome')
            }
            self.samtools_depth = self.add_tool("wgs.samtools_depth")
            self.samtools_depth.set_options(options)
            self.samtools_depth.on("end", self.run_coverage_stat, s)
            self.samtools_depth.run()

    def run_coverage_stat(self, event):
        """
        tool: depth_stat_windows,覆盖度数据处理
        """
        obj = event["bind_object"]
        sample_name = event["data"]
        for f in os.listdir(obj.output_dir):
            depth_file = os.path.join(obj.output_dir, f)
        options = {
            "depth_file": depth_file,
            "step_num": self.option("step_num")
        }
        self.coverage_stat = self.add_tool("wgs.depth_stat_windows")
        self.coverage_stat.set_options(options)
        self.coverage_stat.on("end", self.set_output)
        self.coverage_stat.run()

    def set_output(self, event):
        """
        """
        obj = event["bind_object"]
        insert_path = os.path.join(self.output_dir, "insert")
        depth_path = os.path.join(self.output_dir, "depth")
        coverage_path = os.path.join(self.output_dir, "coverage")
        if not os.path.exists(insert_path):
            os.mkdir(insert_path)
        if not os.path.exists(depth_path):
            os.mkdir(depth_path)
        if not os.path.exists(coverage_path):
            os.mkdir(coverage_path)
        for f in os.listdir(obj.output_dir):
            if f.endswith("insert.xls"):
                if os.path.exists(os.path.join(insert_path, f)):
                    os.remove(os.path.join(insert_path, f))
                os.link(os.path.join(obj.output_dir, f), os.path.join(insert_path, f))
            elif f.endswith("depth.xls"):
                if os.path.exists(os.path.join(depth_path, f)):
                    os.remove(os.path.join(depth_path, f))
                os.link(os.path.join(obj.output_dir, f), os.path.join(depth_path, f))
            elif f.endswith("coverage.xls"):
                if os.path.exists(os.path.join(coverage_path, f)):
                    os.remove(os.path.join(coverage_path, f))
                os.link(os.path.join(obj.output_dir, f), os.path.join(coverage_path, f))
            else:
                if os.path.exists(os.path.join(self.output_dir, f)):
                    os.remove(os.path.join(self.output_dir, f))
                os.link(os.path.join(obj.output_dir, f), os.path.join(self.output_dir, f))
        self.end_times += 1
        if self.start_times == self.end_times:
            time.sleep(1)
            self.get_total_mapped_detail()
            self.end()

    def get_total_mapped_detail(self):
        result_stat = os.path.join(self.output_dir, "result.stat")
        if not os.path.exists(result_stat):
            os.mkdir(result_stat)
        total_mapped = os.path.join(result_stat, "Total.mapped.detail.xls")
        header = "#SampleID\tMapped Ratio(%)\tProperly Mapped(%)\tDuplicate Ratio(%)\tInsert Size\tAverage Depth\t"
        header += "Real Depth\tgenome coverage(1X)\tgenome coverage(5X)\n"
        self.specimen = {}
        for s in self.bam_list.keys():
            name = s.split("-")[0]
            if name not in self.specimen.keys():
                self.specimen[name] = []
            self.specimen[name].append(s)
        with open(total_mapped, "w") as w:
            w.write(header)
            for s in self.specimen.keys():
                map_ratio, prop_map, dup_ratio, insert_size = 0, 0, 0, 0
                average_depth, real_depth, genome_cover1, genome_cover5 = 0, 0, 0, 0
                for s1 in self.specimen[s]:
                    f = os.path.join(self.output_dir, s1 + ".result.stat.xls")
                    with open(f, "r") as r:
                        lines = r.readlines()
                        map_ratio += float(lines[1].strip().split("\t")[1])
                        prop_map += float(lines[2].strip().split("\t")[1])
                        dup_ratio += float(lines[3].strip().split("\t")[1])
                        insert_size += float(lines[4].strip().split("\t")[1])
                        average_depth += float(lines[5].strip().split("\t")[1])
                        real_depth += float(lines[6].strip().split("\t")[1])
                        genome_cover1 += float(lines[-2].strip().split("\t")[1])
                        genome_cover5 += float(lines[-1].strip().split("\t")[1])
                    os.remove(f)
                w.write(s + "\t" + str(map_ratio) + "\t" + str(prop_map) + "\t" + str(dup_ratio) + "\t" + str(insert_size) + "\t"\
                        + str(average_depth) + "\t" + str(real_depth) + "\t" + str(genome_cover1) + "\t" + str(genome_cover5) + "\n")

    def run(self):
        super(MappingStatModule, self).run()
        self.get_bam_list()
        self.run_samtools_stats()
        self.run_samtools_depth()

    def end(self):
        super(MappingStatModule, self).end()

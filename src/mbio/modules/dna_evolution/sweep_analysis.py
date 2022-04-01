# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# last_modify:20180823

import os
import re
import json
from biocluster.module import Module
from biocluster.core.exceptions import OptionError


class SweepAnalysisModule(Module):
    """
    选择性消除
    """
    def __init__(self, work_id):
        super(SweepAnalysisModule, self).__init__(work_id)
        options = [
            {"name": "vcf_file", "type": "infile", "format": "dna_evolution.vcf", "required": True},  # pop.final.vcf/比较分析的结果/vcftools过滤后的结果
            {"name": "group_file", "type": "infile", "format": "dna_evolution.group_table"},  # 分组文件
            {"name": "group_info", "type": "string"},  # 分组方案
            {"name": "diff_group", "type": "string", "required": True},  # 差异分组,逗号分隔,A_vs_B,C_vs_D
            {"name": "recode", "type": "bool", "default": True},
            {"name": "remove_indels", "type": "bool", "default": False},     # --remove-indels
            {"name": "remove_filtered_all", "type": "bool", "default": False},   # --remove-filtered-all
            {"name": "minDP", "type": "int"},           # --minDP 2
            {"name": "maxDP", "type": "int"},           # --maxDP 6
            {"name": "max_missing", "type": "float"},   # --max-missing
            {"name": "min_maf", "type": "float"},       # --maf 0.05
            {"name": "max_maf", "type": "float"},       # --max-maf 0.1
            {"name": "analysis_method", "type": "string", "default": "all", "required": True},  # 分析方法，window_pi、tajima_d、weir_fst_pop
            {"name": "window_size", "type": "int", "default": 2000000},  # --window-pi,窗口大小
            {"name": "window_step", "type": "int", "default": 10000},  # --window-pi-step,窗口步长
            {"name": "vcftools_filter", "type": "bool", "default": True},  # 是否进行vcf过滤，若不过滤，则输入vcf_file是已经过滤过的vcf
            {"name": "update_info", "type": "string"},
            {"name": "main_id", "type": "string"}
        ]
        self.add_option(options)
        self.vcftools_stat, self.draw_stat = [], []

    def check_options(self):
        if not self.option("group_file").is_set and not self.option("group_info"):
            raise OptionError("必需设置分组文件group_file或分组方案group_info")
        if self.option("minDP") and self.option("maxDP"):
            if self.option("minDP") >= self.option("maxDP"):
                raise OptionError("测序深度minDP:{}必须小于minDP:{}".format(self.option("minDP"), self.option("maxDP")))
        if self.option("min_maf") and self.option("max_maf"):
            if self.option("min_maf") >= self.option("max_maf"):
                raise OptionError("次要等位基因频率min_maf:{}必须小于max_maf:{}".format(self.option("min_maf"), self.option("max_maf")))
            if self.option("max_maf") > 1:
                raise OptionError("次要等位基因频率最大值:{}必须小于等于1".format(self.option("max_maf")))
        if self.option("max_missing") < 0 and self.option("max_missing"):
            if self.option("max_missing") < 0 or self.option("max_missing") > 1:
                raise OptionError("缺失率max_missing:{}必须在[0,1]".format(self.option("max_missing")))
        if self.option("window_size") < 0:
            raise OptionError("窗口大小:{}必须大于1".format(self.option("window_size")))
        if self.option("analysis_method") not in ["all", "window_pi", "tajima_d", "weir_fst_pop"]:
            raise OptionError("分析方法:{}只能是all/window_pi/tajima_d/weir_fst_pop,请检查".format(self.option("analysis_method")))

    def get_group_info(self):
        """
        分组文件
        """
        if self.option("group_file").is_set:
            group_info = self.option("group_file").prop["group_info"]
        else:
            group_info = json.loads(self.option("group_info"))
        self.compare_list = self.option("diff_group").split(",")
        self.logger.info(group_info)
        self.logger.info(self.compare_list)
        self.group_list = []
        for diff in self.compare_list:
            group = diff.split("_vs_")[0]
            compare_group = diff.split("_vs_")[1]
            self.group_list.append(group)
            self.group_list.append(compare_group)
        self.group_list = list(set(self.group_list))
        self.group_paths = {}
        for group in self.group_list:
            group_path = self.work_dir + "/" + group + ".list"
            self.group_paths[group] = group_path
            with open(group_path, "w") as w:
                w.write("\n".join(group_info[group]))

    def run_vcftools_filter(self):
        """
        对vcf文件进行过滤
        """
        options = {
            "vcf_path": self.option("vcf_file").prop["path"],
            "recode": self.option("recode"),
            "remove_indels": self.option("remove_indels"),
            "remove_filtered_all": self.option("remove_filtered_all")
        }
        if self.option("minDP"):
            options["minDP"] = self.option("minDP")
        if self.option("maxDP"):
            options["maxDP"] = self.option("maxDP")
        if self.option("max_missing"):
            options["max_missing"] = self.option("max_missing")
        if self.option("min_maf"):
            options["min_maf"] = self.option("min_maf")
        if self.option("max_maf"):
            options["max_maf"] = self.option("max_maf")
        self.vcftools_filter = self.add_tool("dna_evolution.vcftools_filter")
        self.vcftools_filter.set_options(options)
        self.vcftools_filter.on("end", self.check_vcf)
        self.vcftools_filter.on("end", self.set_output, "vcf_filter")
        self.vcftools_filter.run()

    def run_vcftools_stat_tajimad(self):
        for group in self.group_list:
            options = {
                "group_list": self.group_paths[group],
                "analysis_method": "tajima_d",
                "window_size": self.option("window_size"),
                "window_step": self.option("window_step"),
            }
            if self.option("vcftools_filter"):
                options["vcf_file"] = self.vcftools_filter.option("filter_vcf").prop["path"]
            else:
                options["vcf_file"] = self.option("vcf_file").prop["path"]
            self.vcftools_tajimad = self.add_tool("dna_evolution.vcftools_stat")
            self.vcftools_tajimad.set_options(options)
            self.vcftools_tajimad.on("end", self.set_output, "tajimad_" + group)
            self.vcftools_stat.append(self.vcftools_tajimad)
        self.set_on_rely("vcftools_stat")

    def run_vcftools_stat_window_pi(self):
        for group in self.group_list:
            options = {
                "group_list": self.group_paths[group],
                "analysis_method": "window_pi",
                "window_size": self.option("window_size"),
                "window_step": self.option("window_step"),
            }
            if self.option("vcftools_filter"):
                options["vcf_file"] = self.vcftools_filter.option("filter_vcf").prop["path"]
            else:
                options["vcf_file"] = self.option("vcf_file").prop["path"]
            self.vcftools_window_pi = self.add_tool("dna_evolution.vcftools_stat")
            self.vcftools_window_pi.set_options(options)
            self.vcftools_window_pi.on("end", self.set_output, "window_pi_" + group)
            self.vcftools_stat.append(self.vcftools_window_pi)
        self.set_on_rely("vcftools_stat")

    def run_vcftools_stat_weir_fst_pop(self):
        self.logger.info(self.compare_list)
        for cg in self.compare_list:
            group = cg.split("_vs_")[0]
            compare_group = cg.split("_vs_")[1]
            options = {
                "group_list": self.group_paths[group],
                "compare_group": self.group_paths[compare_group],
                "analysis_method": "weir_fst_pop",
                "window_size": self.option("window_size"),
                "window_step": self.option("window_step"),
            }
            if self.option("vcftools_filter"):
                options["vcf_file"] = self.vcftools_filter.option("filter_vcf").prop["path"]
            else:
                options["vcf_file"] = self.option("vcf_file").prop["path"]
            self.vcftools_weir_fst_pop = self.add_tool("dna_evolution.vcftools_stat")
            self.vcftools_weir_fst_pop.set_options(options)
            self.vcftools_weir_fst_pop.on("end", self.set_output, "weir_fst_pop_" + group)
            self.vcftools_stat.append(self.vcftools_weir_fst_pop)
        self.set_on_rely("vcftools_stat")

    def set_on_rely(self, type):
        if type == "vcftools_stat":
            if len(self.vcftools_stat) == 2 * len(self.group_list) + len(self.compare_list):
                self.on_rely(self.vcftools_stat, self.run_draw_data)
                self.logger.info("设定第一次依赖成功")
                for tool in self.vcftools_stat:
                    tool.run()
        else:
            if len(self.draw_stat) == len(self.group_list) + len(self.compare_list):
                self.on_rely(self.draw_stat, self.set_db)
                self.logger.info("设定第二层依赖成功")
                for tool in self.draw_stat:
                    tool.run()

    def run_draw_data(self):
        self.run_draw_pi_tajima()
        self.run_fst_manhattan()

    def run_draw_pi_tajima(self):
        self.draw_pi_tajima_stats = []
        for group in self.group_list:
            options = {
                "tajima_file": self.output_dir + "/" + group + ".Tajima.D",
                "window_pi_file": self.output_dir + "/" + group + ".windowed.pi"
            }
            self.draw_pi_tajima = self.add_tool("dna_evolution.draw_pi_tajima")
            self.draw_pi_tajima.set_options(options)
            # self.draw_pi_tajima.on("end", self.set_output, "draw_pi_tajima_" + group)
            self.draw_stat.append(self.draw_pi_tajima)
        self.set_on_rely("draw_pi_tajima")

    def run_fst_manhattan(self):
        self.draw_fst_manhattan_stats = []
        for cg in self.compare_list:
            group = cg.split("_vs_")[0]
            compare_group = cg.split("_vs_")[1]
            options = {
                "tajima": self.output_dir + "/" + group + ".Tajima.D",
                "window_pi": self.output_dir + "/" + group + ".windowed.pi",
                "compare_tajima": self.output_dir + "/" + compare_group + ".Tajima.D",
                "compare_window_pi": self.output_dir + "/" + compare_group + ".windowed.pi",
                "window_weir_fst": self.output_dir + "/" + group + "-" + compare_group + ".windowed.weir.fst"
            }
            self.fst_manhattan = self.add_tool("dna_evolution.fst_manhattan")
            self.fst_manhattan.set_options(options)
            self.fst_manhattan.on("end", self.set_output, "fst_manhattan_" + group)
            self.draw_stat.append(self.fst_manhattan)
        self.set_on_rely("fst_manhattan")

    def get_vcf_variant_num(self):
        chr_list = []
        window_step = self.option("window_step")
        if self.option("vcftools_filter"):
            vcf_path = self.vcftools_filter.option("filter_vcf").prop["path"]
        else:
            vcf_path = self.option("vcf_file").prop["path"]
        vcf_variant_file = os.path.join(self.work_dir, "vcf.variant.txt")
        with open(vcf_path, "r") as r, open(vcf_variant_file, "w") as w:
            w.write("#Chr\tPos\tVariant Num\n")
            lines = r.readlines()
            for line in lines:
                if line.startswith("#"):
                    continue
                item = line.strip().split("\t")
                chr = item[0]
                pos = int(item[1])
                if item[0] not in chr_list:
                    chr_list.append(chr)
                    start = pos / window_step * window_step + 1
                    min = start
                    variant_num = 0
                max = min + window_step
                if pos >= min and pos < max:
                    variant_num += 1
                else:
                    w.write(chr + "\t" + str(min) + "\t" + str(variant_num) + "\n")
                    more_num = (pos - max) / window_step
                    if more_num > 0:
                        # for i in range(more_num):  # 当没有这个范围内的pos的时候写入0
                        #     min = max + window_step * i
                        #     w.write(chr + "\t" + str(min) + "\t" + str(0) + "\n")
                        min = max + window_step * more_num
                    else:
                        min = max
                    variant_num = 1
            w.write(chr + "\t" + str(min) + "\t" + str(variant_num) + "\n")

    def set_output(self, event):
        obj = event["bind_object"]
        if event["data"].startswith("fst_manhattan_"):
            for f in os.listdir(obj.output_dir):
                if f.endswith("pi_tajimaD_fst.detail"):
                    self.link_file(os.path.join(obj.output_dir, f), os.path.join(self.output_dir, f))
                if re.search("png|pdf", f):
                    self.link_file(os.path.join(obj.output_dir, f), os.path.join(self.output_dir, f))

        else:
            for f in os.listdir(obj.output_dir):
                if re.search(r"pdf|png|select", f):
                    continue
                self.link_file(os.path.join(obj.output_dir, f), os.path.join(self.output_dir, f))

    def link_file(self, old, new):
        if os.path.exists(new):
            os.remove(new)
        os.link(old, new)

    def set_db(self):
        if self.option("main_id"):
            self.logger.info("将选择性消除结果导入mongo表")
            sweep_api = self.api.api("dna_evolution.sweep_analysis")
            sweep_id = self.option("main_id")
            window_step = self.option("window_step")
            variant_num_path = os.path.join(self.work_dir, "vcf.variant.txt")
            sweep_dir = self._sheet.output + "/"
            for diff_group_name in self.option("diff_group").split(","):
                pop1 = diff_group_name.split("_vs_")[0]
                pop2 = diff_group_name.split("_vs_")[1]
                pi_tajimad_fst_path = os.path.join(self.output_dir, pop1 + "-" + pop2 + ".pi_tajimaD_fst.detail")
                if os.path.exists(pi_tajimad_fst_path):
                    sweep_api.add_sg_sweep_detail(sweep_id, pi_tajimad_fst_path, variant_num_path, diff_group_name, window_step)
                else:
                    self.bind_object.logger.info("差异分组：{}没有pi_tajimaD_fst结果".format(diff_group_name))
                png_path = os.path.join(sweep_dir, pop1 + "-" + pop2 + ".pop1.manhattan.png")
                pdf_path = os.path.join(sweep_dir, pop1 + "-" + pop2 + ".pop1.manhattan.pdf")
                sweep_api.add_sg_manhattan_path(sweep_id, diff_group_name, pop1, png_path, pdf_path)
                png_path = os.path.join(sweep_dir, pop1 + "-" + pop2 + ".pop2.manhattan.png")
                pdf_path = os.path.join(sweep_dir, pop1 + "-" + pop2 + ".pop2.manhattan.pdf")
                sweep_api.add_sg_manhattan_path(sweep_id, diff_group_name, pop2, png_path, pdf_path)
            vcf_path = os.path.join(sweep_dir, "pop.recode.vcf")
            sweep_api.update_sweep_path(sweep_id, sweep_dir, vcf_path)
        self.end()

    def check_vcf(self):
        """
        检查vcf文件是否是空的
        """
        if self.option("vcftools_filter"):
            vcf_path = self.vcftools_filter.option("filter_vcf").prop["path"]
        else:
            vcf_path = self.option("vcf_file").prop["path"]
        with open(vcf_path, "r") as f:
            num = 0
            lines = f.readline().strip()
            while lines:
                if not re.match('#', lines):
                    num = 1     # 筛查只有#开头的内容,不然
                    break       # 存在非#开头行的话，break
                lines = f.readline().strip()
            if num == 0:
                self.logger.info("vcf文件:%s是空的" % vcf_path)
                self.end()
        self.run_vcftools_stat_tajimad()
        self.run_vcftools_stat_window_pi()
        self.run_vcftools_stat_weir_fst_pop()
        self.get_vcf_variant_num()

    def run(self):
        super(SweepAnalysisModule, self).run()
        self.get_group_info()
        if self.option("vcftools_filter"):
            self.run_vcftools_filter()
        else:
            self.check_vcf()
            # self.run_vcftools_stat_tajimad()
            # self.run_vcftools_stat_window_pi()
            # self.run_vcftools_stat_weir_fst_pop()
            # self.get_vcf_variant_num()

    def end(self):
        self.add_upload_dir(self.output_dir)
        super(SweepAnalysisModule, self).end()

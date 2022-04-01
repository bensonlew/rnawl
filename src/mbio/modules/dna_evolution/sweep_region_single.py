# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# last_modify:20180919

import os
import re
from biocluster.module import Module
from biocluster.core.exceptions import OptionError


class SweepRegionSingleModule(Module):
    """
    受选择区域筛选
    """
    def __init__(self, work_id):
        super(SweepRegionSingleModule, self).__init__(work_id)
        options = [
            {"name": "vcf_file", "type": "infile", "format": "dna_gmap.vcf", "required": True},  # pop.recode.vcf或者比较分析的结果
            {"name": "pop_summary", "type": "infile", "format": "dna_gmap.vcf", "required": True},  # pop.summary
            {"name": "tajima", "type": "infile", "format": "dna_evolution.tajima_d", "required": True},  # tajima_d
            {"name": "window_pi", "type": "infile", "format": "dna_evolution.window_pi", "required": True},  # window_pi
            {"name": "compare_tajima", "type": "infile", "format": "dna_evolution.tajima_d"},  # tajima_d
            {"name": "compare_window_pi", "type": "infile", "format": "dna_evolution.window_pi"},  # window_pi
            {"name": "window_weir_fst", "type": "infile", "format": "dna_evolution.window_weir_fst"},  # window_weir_fst
            {"name": "p_value", "type": "float", "default": 0.05},  # manhattan的p_value筛选
            {"name": "ld_window", "type": "int", "default": 99999},  # --ld-window
            {"name": "ld_window_kb", "type": "int", "default": 500},  # --ld-window-kb, 固定大小
            {"name": "ld_window_r2", "type": "float", "default": 0.8},  # --ld-window-r2, R2
            {"name": "chr_set", "type": "int", "default": 35},   # --chr-set 35
            {"name": "pop_type", "type": "string", "default": "diff_pop"},  # 群体类型,diff_pop/pop
        ]
        self.add_option(options)

    def check_options(self):
        if self.option("p_value") <= 0 or self.option("p_value") >= 1:
            raise OptionError("p_value:%s的范围是0-1,请检查" % self.option("p_value"))
        if self.option("ld_window_kb") <= 0:
            raise OptionError("固定大小：{}必须大于0".format(self.option("ld_window_kb")))
        if self.option("ld_window_r2") > 1 or self.option("ld_window_r2") < 0:
            raise OptionError("R2的大小范围是0-1，而不是:{}".format(self.option("ld_window_r2")))
        if self.option("chr_set") <= 0:
            raise OptionError("染色体:{}必须大于0".format(self.option("chr_set")))
        if self.option("pop_type") not in ["diff_pop", "pop"]:
            raise OptionError("pop_type:{}只能是diff_pop/pop,请检查")
        if self.option("pop_type") == "diff_pop":
            if not self.option("compare_window_pi").is_set:
                raise OptionError("pop_type为diff_pop的时候必须设置compare_window_pi文件")
            if not self.option("compare_tajima").is_set:
                raise OptionError("pop_type为diff_pop的时候必须设置compare_tajima文件")
            if not self.option("window_weir_fst").is_set:
                raise OptionError("pop_type为diff_pop的时候必须设置window_weir_fst文件")

    def run_pi_tajimad(self):
        options = {
            "tajima_file": self.option("tajima"),
            "window_pi_file": self.option("window_pi"),
            "p_value": self.option("p_value")
        }
        self.pi_tajimad = self.add_tool("dna_evolution.draw_pi_tajima")
        self.pi_tajimad.set_options(options)
        self.pi_tajimad.on("end", self.run_region_snp)
        self.pi_tajimad.run()

    def run_fst_pi(self):
        options = {
            "window_pi": self.option("window_pi"),
            "compare_window_pi": self.option("compare_window_pi"),
            "window_weir_fst": self.option("window_weir_fst"),
            "p_value": self.option("p_value"),
            "analysis_type": "fst_pi"
        }
        self.fst_pi = self.add_tool("dna_evolution.fst_manhattan")
        self.fst_pi.set_options(options)
        self.fst_pi.on("end", self.run_region_snp)
        self.fst_pi.run()

    def run_region_snp(self):
        options = {"vcf_file": self.option("vcf_file")}
        if self.option("pop_type") == "pop":
            if self.pi_tajimad.option("select_file").prop["is_null"]:
                self.end()
            else:
                options["select_file"] = self.pi_tajimad.option("select_file").prop["path"]
                self.region_snp = self.add_tool("dna_evolution.region_snp")
                self.region_snp.set_options(options)
                self.region_snp.on("end", self.run_vcftools_plink)
                self.region_snp.run()
        else:
            if self.fst_pi.option("select_file").prop["is_null"]:
                self.end()
            else:
                options["select_file"] = self.fst_pi.option("select_file").prop["path"]
                self.region_snp = self.add_tool("dna_evolution.region_snp")
                self.region_snp.set_options(options)
                self.region_snp.on("end", self.run_vcftools_plink)
                self.region_snp.run()

    def run_vcftools_plink(self):
        options = {
            "r2": True,
            "flag": False,
            "allow_extra_chr": True,
            "chr_set": self.option("chr_set"),
            "recode_vcf_path": self.option("vcf_file"),
            "ldsnp_list": self.region_snp.option("region_snp")
        }
        self.vcftools_plink = self.add_tool("dna_evolution.vcftools_plink")
        self.vcftools_plink.set_options(options)
        self.vcftools_plink.on("end", self.run_ld_region)
        self.vcftools_plink.run()

    def run_ld_region(self):
        for f in os.listdir(self.vcftools_plink.output_dir):
            if f.endswith("select.snp.ld"):
                snp_ld = os.path.join(self.vcftools_plink.output_dir, f)
        if not snp_ld:
            self.set_error("没有在vcftools的output里找到1.pi_tajimaD_fst.select.snp.ld结尾的文件")
        options = {
            "snp_ld": snp_ld
        }
        self.ld_region = self.add_tool("dna_evolution.ld_region")
        self.ld_region.set_options(options)
        self.ld_region.on("end", self.run_sweep_region_stat)
        self.ld_region.on("end", self.set_output, "ld_region")
        self.ld_region.run()

    def run_sweep_region_stat(self):
        for f in os.listdir(self.ld_region.output_dir):
            region_file = os.path.join(self.ld_region.output_dir, f)
        if not region_file:
            self.set_error("没有在ld_region的output里找到文件")
        options = {
            "vcf_file": self.option("vcf_file"),
            "pop_summary": self.option("pop_summary"),
            "tajima": self.option("tajima"),
            "window_pi": self.option("window_pi"),
            "region_file": region_file,
            "pop_type": self.option("pop_type")
        }
        if self.option("pop_type") == "diff_pop":
            options["window_weir_fst"] = self.option("window_weir_fst")
            options["compare_tajima"] = self.option("compare_tajima")
            options["compare_window_pi"] = self.option("compare_window_pi")
        self.sweep_region_stat = self.add_tool("dna_evolution.sweep_region_stat")
        self.sweep_region_stat.set_options(options)
        self.sweep_region_stat.on("end", self.set_output, "sweep_region_stat")
        self.sweep_region_stat.run()

    def set_output(self, event):
        obj = event["bind_object"]
        for f in os.listdir(obj.output_dir):
            new = os.path.join(self.output_dir, f)
            old = os.path.join(obj.output_dir, f)
            if os.path.exists(new):
                os.remove(new)
            os.link(old, new)
        if event["data"] == "sweep_region_stat":
            self.end()

    def run(self):
        super(SweepRegionSingleModule, self).run()
        if self.option("pop_type") == "pop":
            self.run_pi_tajimad()
        else:
            self.run_fst_pi()

    def end(self):
        super(SweepRegionSingleModule, self).end()

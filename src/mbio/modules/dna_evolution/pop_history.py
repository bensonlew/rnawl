# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# last_modify:20180921

import os
import json
from biocluster.module import Module
from biocluster.core.exceptions import OptionError


class PopHistoryModule(Module):
    """
    种群历史
    """
    def __init__(self, work_id):
        super(PopHistoryModule, self).__init__(work_id)
        options = [
            {"name": "vcf_file", "type": "infile", "format": "dna_gmap.vcf", "required": True},  # pop.final.vcf/比较分析的结果/vcftools过滤后的结果
            {"name": "group_file", "type": "infile", "format": "dna_evolution.group_table"},  # 分组文件
            {"name": "group_info", "type": "string"},  # 分组方案
            {"name": "recode", "type": "bool", "default": True},
            {"name": "remove_indels", "type": "bool", "default": False},  # --remove-indels
            {"name": "remove_filtered_all", "type": "bool", "default": False},  # --remove-filtered-all
            {"name": "minDP", "type": "int"},  # --minDP 1
            {"name": "maxDP", "type": "int"},  # --maxDP 30
            {"name": "max_missing", "type": "float"},  # --max-missing 0.3
            {"name": "min_maf", "type": "float"},  # --maf 0.05
            {"name": "max_maf", "type": "float"},  # --max-maf 1
            {"name": "iteration_num", "type": "int", "default": 25},  # -N,最大迭代次数
            {"name": "coalescent_time", "type": "float", "default": 25},  # -t,最大2NO合并时间
            {"name": "theta_ratio", "type": "float", "default": 5},  # -r,初始的theta/rho比率
            {"name": "parameters_pattern", "type": "string", "default": "4+25*2+4+6"},  # -p,参数模式
            {"name": "analysis_method", "type": "string", "default": "psmc"},  # 分析方法，psmc,smc
            {"name": "vcftools_filter", "type": "bool", "default": True},  # 是否进行vcf过滤，若不过滤，则输入vcf_file是已经过滤过的vcf
            {"name": "update_info", "type": "string"},
            {"name": "main_id", "type": "string"},
        ]
        self.add_option(options)
        self.all_modules = []

    def check_options(self):
        if not self.option("group_file").is_set and not self.option("group_info"):
            raise OptionError("必需设置分组文件group_file或分组方案group_info")
        if self.option("minDP") and self.option("maxDP"):
            if self.option("minDP") >= self.option("maxDP"):
                raise OptionError("测序深度minDP:{}必须小于maxDP:{}".format(self.option("minDP"), self.option("maxDP")))
        if self.option("min_maf") and self.option("max_maf"):
            if self.option("min_maf") >= self.option("max_maf"):
                raise OptionError("次要等位基因频率min_maf:{}必须小于max_maf:{}".format(self.option("min_maf"), self.option("max_maf")))
        if self.option("max_maf"):
            if self.option("max_maf") > 1:
                raise OptionError("次要等位基因频率最大值:{}必须小于等于1".format(self.option("max_maf")))
        if self.option("max_missing"):
            if self.option("max_missing") < 0 or self.option("max_missing") > 1:
                raise OptionError("缺失率max_missing:{}必须在[0,1]".format(self.option("max_missing")))
        if self.option("analysis_method") not in ["psmc", "smc"]:
            raise OptionError("分析方法:{}只能是psmc/smc".format(self.option("analysis_method")))

    def get_group_info(self):
        """
        分组文件
        """
        if self.option("group_file").is_set:
            self.group_file = self.option("group_file").prop["path"]
        else:
            self.group_file = os.path.join(self.work_dir, "group.list")
            group_info = json.loads(self.option("group_info"))
            with open(self.group_file, "w") as w:
                for group in group_info.keys():
                    for sample in group_info[group]:
                        w.write(sample + "\t" + group + "\n")
                        w.write(sample + "\t" + "all" + "\n")

    def run_vcftools_filter(self):
        """
        对vcf文件进行过滤
        """
        options = {
            "vcf_path": self.option("vcf_file"),
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
        self.vcftools_filter.set_options(options)
        self.vcftools_filter.run()

    def run_vcf2psmc(self):
        """
        得到每个组的fa文件
        """
        options = {
            "group_list": self.group_file
        }
        if self.option("vcftools_filter"):
            options["vcf_file"] = self.vcftools_filter.option("filter_vcf").prop["path"]
        else:
            options["vcf_file"] = self.option("vcf_file").prop["path"]
        self.vcf2psmc = self.add_tool("dna_evolution.vcf2psmc")
        self.vcf2psmc.set_options(options)
        self.vcf2psmc.on("end", self.run_psmc)
        self.vcf2psmc.run()

    def run_psmc(self):
        """
        psmc
        """
        options = {
            "psmcfa_list": self.vcf2psmc.option("psmcfa_list"),
            "iteration_num": self.option("iteration_num"),
            "coalescent_time": self.option("coalescent_time"),
            "theta_ratio": self.option("theta_ratio"),
            "parameters_pattern": self.option("parameters_pattern")
        }
        self.psmc = self.add_tool("dna_evolution.psmc")
        self.psmc.set_options(options)
        self.psmc.on("end", self.run_psmc_stat)
        self.psmc.run()

    def run_psmc_stat(self):
        """
        psmc
        """
        options = {
            "psmc_list": self.psmc.option("psmc_list"),
        }
        self.psmc_stat = self.add_tool("dna_evolution.psmc_stat")
        self.psmc_stat.set_options(options)
        self.psmc_stat.on("end", self.set_output, "psmc_stat")
        self.psmc_stat.run()

    def run_vcf2smc(self):
        """
        vcf2smc
        """
        options = {
            "group_list": self.group_file
        }
        if self.option("vcftools_filter"):
            options["vcf_file"] = self.vcftools_filter.option("filter_vcf").prop["path"]
        else:
            options["vcf_file"] = self.option("vcf_file").prop["path"]
        self.vcf2smc = self.add_tool("dna_evolution.vcf2smc")
        self.vcf2smc.set_options(options)
        self.vcf2smc.on("end", self.run_smc)
        self.vcf2smc.run()

    def run_smc(self):
        """
        smc
        """
        self.all_smc = []
        for group in os.listdir(self.vcf2smc.output_dir):
            options = {
                "smc_dir": os.path.join(self.vcf2smc.output_dir, group)
            }
            self.smc = self.add_tool("dna_evolution.smc")
            self.smc.set_options(options)
            self.smc.on("end", self.get_smc_list, "smc_list")
            self.all_smc.append(self.smc)
        if len(self.all_smc) == 1:
            self.all_smc[0].on("end", self.run_smc_stat)
            self.all_smc[0].run()
        else:
            self.on_rely(self.all_smc, self.run_smc_stat)
            for smc in self.all_smc:
                smc.run()

    def run_smc_stat(self):
        """
        smc_stat
        """
        options = {
            "smc_dir": os.path.join(self.work_dir, "smc_result")
        }
        self.smc_stat = self.add_tool("dna_evolution.smc_stat")
        self.smc_stat.set_options(options)
        self.smc_stat.on("end", self.set_output, "smc_stat")
        self.smc_stat.run()

    def get_smc_list(self, event):
        """
        将smc的结果link到smc_result文件夹里
        """
        obj = event["bind_object"]
        smc_result = os.path.join(self.work_dir, "smc_result")
        if not os.path.exists(smc_result):
            os.mkdir(smc_result)
        for f in os.listdir(obj.output_dir):
            if f.endswith("smc.result"):
                new = os.path.join(smc_result, f.split(".smc.result")[0] + ".psmc.result")
                old = os.path.join(obj.output_dir, f)
                if os.path.exists(new):
                    os.remove(new)
                os.link(old, new)

    def set_output(self, event):
        obj = event["bind_object"]
        for f in os.listdir(obj.output_dir):
            new = os.path.join(self.output_dir, f)
            old = os.path.join(obj.output_dir, f)
            if os.path.exists(new):
                os.remove(new)
            os.link(old, new)
        self.set_db()

    def set_db(self):
        if self.option("main_id"):
            self.logger.info("将种群历史结果导入mongo表")
            history_api = self.api.api("dna_evolution.pop_history")
            history_id = self.option("main_id")
            if self.option("analysis_method") == "psmc":
                stat_path = os.path.join(self.output_dir, "all_psmc_stat.xls")
            else:
                stat_path = os.path.join(self.output_dir, "all_smc_stat.xls")
            history_api.add_sg_history_detail(history_id, stat_path)
            curve_id = history_api.add_sg_curve(history_id)
            history_api.add_sg_curve_detail(curve_id, stat_path)
        self.end()

    def run(self):
        super(PopHistoryModule, self).run()
        self.get_group_info()
        if self.option("vcftools_filter"):
            self.vcftools_filter = self.add_tool("dna_evolution.vcftools_filter")
            if self.option("analysis_method") == "psmc":
                self.vcftools_filter.on("end", self.run_vcf2psmc)
            else:
                self.vcftools_filter.on("end", self.run_vcf2smc)
            self.run_vcftools_filter()
        else:
            if self.option("analysis_method") == "psmc":
                self.run_vcf2psmc()
            else:
                self.run_vcf2smc()

    def end(self):
        super(PopHistoryModule, self).end()

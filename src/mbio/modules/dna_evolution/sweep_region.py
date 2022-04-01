# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# last_modify:20180919

import os
from biocluster.module import Module
from biocluster.core.exceptions import OptionError


class SweepRegionModule(Module):
    """
    受选择区域筛选
    """
    def __init__(self, work_id):
        super(SweepRegionModule, self).__init__(work_id)
        options = [
            {"name": "vcf_file", "type": "infile", "format": "dna_gmap.vcf", "required": True},  # 选择性消除已经进行vcf过滤后的vcf文件
            {"name": "pop_summary", "type": "infile", "format": "dna_gmap.vcf", "required": True},  # pop.summary
            {"name": "sweep_dir", "type": "infile", "format": "dna_evolution.sweep_dir", "required": True},  # 选择性消除的结果
            {"name": "diff_group", "type": "string", "required": True},  # 差异分组,逗号分隔
            {"name": "p_value", "type": "float", "default": 0.05},  # manhattan的p_value筛选
            {"name": "ld_window", "type": "int", "default": 99999},  # --ld-window
            {"name": "ld_window_kb", "type": "int", "default": 500},  # --ld-window-kb, 固定大小
            {"name": "ld_window_r2", "type": "float", "default": 0.8},  # --ld-window-r2, R2
            {"name": "chr_set", "type": "int", "default": 35},   # --chr-set 35
            {"name": "update_info", "type": "string"},
            {"name": "main_id", "type": "string"}
        ]
        self.add_option(options)
        self.all_modules = []

    def check_options(self):
        if self.option("p_value") <= 0 or self.option("p_value") >= 1:
            raise OptionError("p_value:%s的范围是0-1,请检查" % self.option("p_value"))
        if self.option("ld_window_kb") <= 0:
            raise OptionError("固定大小：{}必须大于0".format(self.option("ld_window_kb")))
        if self.option("ld_window_r2") > 1 or self.option("ld_window_r2") < 0:
            raise OptionError("R2的大小范围是0-1，而不是:{}".format(self.option("ld_window_r2")))
        if self.option("chr_set") <= 0:
            raise OptionError("染色体:{}必须大于0".format(self.option("chr_set")))

    def get_diff_group(self):
        """
        得到差异分组对应的文件
        """
        self.diff_group, self.group_info = {}, {}
        self.group_list = []
        for diff in self.option("diff_group").split(","):
            group = diff.split("_vs_")[0]
            compare_group = diff.split("_vs_")[1]
            tajima = os.path.join(self.option("sweep_dir").prop["path"], group + ".Tajima.D")
            compare_tajima = os.path.join(self.option("sweep_dir").prop["path"], compare_group + ".Tajima.D")
            window_pi = os.path.join(self.option("sweep_dir").prop["path"], group + ".windowed.pi")
            compare_window_pi = os.path.join(self.option("sweep_dir").prop["path"], compare_group + ".windowed.pi")
            window_weir_fst = os.path.join(self.option("sweep_dir").prop["path"], group + "-" + compare_group + ".windowed.weir.fst")
            self.diff_group[diff] = [tajima, compare_tajima, window_pi, compare_window_pi, window_weir_fst]
            if group not in self.group_list:
                self.group_list.append(group)
                self.group_info[group] = [tajima, window_pi]
            if compare_group not in self.group_list:
                self.group_list.append(compare_group)
                self.group_info[compare_group] = [tajima, window_pi]
        self.logger.info(self.diff_group)
        self.logger.info(self.group_info)

    def run_sweep_region_diff(self):
        """
        每个差异分群
        """
        for diff in self.option("diff_group").split(","):
            options = {
                "vcf_file": self.option("vcf_file"),
                "pop_summary": self.option("pop_summary"),
                "p_value": self.option("p_value"),
                "ld_window": self.option("ld_window"),
                "ld_window_kb": self.option("ld_window_kb"),
                "ld_window_r2": self.option("ld_window_r2"),
                "chr_set": self.option("chr_set"),
                "tajima": self.diff_group[diff][0],
                "window_pi": self.diff_group[diff][2],
                "compare_tajima": self.diff_group[diff][1],
                "compare_window_pi": self.diff_group[diff][3],
                "window_weir_fst": self.diff_group[diff][4],
                "pop_type": "diff_pop"
            }
            self.sweep_region_single = self.add_module("dna_evolution.sweep_region_single")
            self.sweep_region_single.set_options(options)
            self.sweep_region_single.on("end", self.set_output, "sweep_diff_pop")
            self.all_modules.append(self.sweep_region_single)

    def run_sweep_region_single_pop(self):
        """
        每个分群
        """
        for diff in self.group_list:
            options = {
                "vcf_file": self.option("vcf_file"),
                "pop_summary": self.option("pop_summary"),
                "p_value": self.option("p_value"),
                "ld_window": self.option("ld_window"),
                "ld_window_kb": self.option("ld_window_kb"),
                "ld_window_r2": self.option("ld_window_r2"),
                "chr_set": self.option("chr_set"),
                "tajima": self.group_info[diff][0],
                "window_pi": self.group_info[diff][1],
                "pop_type": "pop"
            }
            self.sweep_region_single = self.add_module("dna_evolution.sweep_region_single")
            self.sweep_region_single.set_options(options)
            self.sweep_region_single.on("end", self.set_output, "sweep_single_pop")
            self.all_modules.append(self.sweep_region_single)

    def set_output(self, event):
        obj = event["bind_object"]
        for f in os.listdir(obj.output_dir):
            new = os.path.join(self.output_dir, f)
            old = os.path.join(obj.output_dir, f)
            if os.path.exists(new):
                os.remove(new)
            os.link(old, new)

    def set_db(self):
        if self.option("main_id"):
            self.logger.info("将选择性消除结果导入mongo表")
            region_api = self.api.api("dna_evolution.sweep_region")
            sweep_region_id = self.option("main_id")
            diff_group = self.option("diff_group").split(",")
            for diff_group_name in diff_group:
                name = diff_group_name.split("_vs_")[0] + "-" + diff_group_name.split("_vs_")[1]
                select_region_stat = os.path.join(self.output_dir, name + ".fst_pi.detail.select.snp.ld.snp.stat")
                if os.path.exists(select_region_stat):
                    region_api.add_sg_sweep_region_detail(sweep_region_id, diff_group_name, select_region_stat)
                else:
                    self.logger.info("区域筛选结果为空，请根据具体情况重新调整参数")
            sweep_region_dir = self._sheet.output
            region_api.update_diff_group(sweep_region_id, diff_group, sweep_region_dir)
        self.end()

    def run(self):
        super(SweepRegionModule, self).run()
        self.get_diff_group()
        self.run_sweep_region_diff()
        self.run_sweep_region_single_pop()
        if len(self.all_modules) == 1:
            self.all_modules[0].on("end", self.set_db)
            self.all_modules[0].run()
        else:
            self.on_rely(self.all_modules, self.set_db)
            for module in self.all_modules:
                module.run()

    def end(self):
        self.add_upload_dir(self.output_dir)
        super(SweepRegionModule, self).end()

# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# modified 2018.08.15

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os


class VcftoolsStatAgent(Agent):
    """
    工具：vcftools
    对vcf进行window-pi、TajimaD、weir-fst-pop统计
    """
    def __init__(self, parent):
        super(VcftoolsStatAgent, self).__init__(parent)
        options = [
            {"name": "vcf_file", "type": "infile", "format": "dna_gmap.vcf", "required": True},  # pop.recode.vcf或者比较分析的结果
            {"name": "group_list", "type": "string", "required": True},  # 分组文件,输出文件以分组文件名称.分割命名
            {"name": "compare_group", "type": "string"},  # 方法为weir_fst_pop时的分组文件
            {"name": "analysis_method", "type": "string", "required": True},  # 分析方法，window_pi、tajima_d、weir_fst_pop
            {"name": "window_size", "type": "int", "default": 2000000},  # --window-pi,窗口大小
            {"name": "window_step", "type": "int", "default": 10000},  # --window-pi-step,窗口步长
            {"name": "stat_file", "type": "outfile", "format": "dna_evolution.window_pi,dna_evolution.tajima_d,dna_evolution.window_weir_fst"},  # vcftools进行统计后的输出
        ]
        self.add_option(options)

    def check_options(self):
        if self.option("analysis_method") not in ["window_pi", "tajima_d", "weir_fst_pop"]:
            raise OptionError("参数analysis_method:{}必须为window_pi/tajima_d/weir_fst_pop")
        if self.option("analysis_method") == "weir_fst_pop":
            if not self.option("compare_group"):
                raise OptionError("缺少参数compare_group，请检查")
        with open(self.option("group_list"), "r") as f:
            lines = f.readlines()

    def set_resource(self):
        self._cpu = 2
        self._memory = "15G"

    def end(self):
        super(VcftoolsStatAgent, self).end()


class VcftoolsStatTool(Tool):
    def __init__(self, config):
        super(VcftoolsStatTool, self).__init__(config)
        self.set_environ(LD_LIBRARY_PATH=self.config.SOFTWARE_DIR + '/gcc/5.4.0/lib64')
        self.vcftools = "bioinfo/dna_evolution/vcftools"

    def run_vcftools_window_pi(self):
        """
        window-pi方法
        """
        self.stat_path = self.output_dir + "/" + os.path.basename(self.option("group_list")).split(".")[0]
        cmd = "{} --vcf {} --keep {}".format(self.vcftools, self.option("vcf_file").prop["path"], self.option("group_list"))
        cmd += " --window-pi {} --window-pi-step {} --out {}".format(self.option("window_size"), self.option("window_step"), self.stat_path)
        command = self.add_command("vcftools_window_pi", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("vcftools_window_pi运行成功")
        else:
            self.set_error("vcftools_window_pi运行失败")
        self.option("stat_file", self.stat_path + ".windowed.pi")

    def run_vcftools_tajima_d(self):
        """
        TajimaD方法
        """
        self.stat_path = self.output_dir + "/" + os.path.basename(self.option("group_list")).split(".")[0]
        cmd = "{} --vcf {} --keep {}".format(self.vcftools, self.option("vcf_file").prop["path"], self.option("group_list"))
        cmd += " --TajimaD {} --out {}".format(self.option("window_step"), self.stat_path)
        command = self.add_command("vcftools_tajimad", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("vcftools_tajimad运行成功")
        else:
            self.set_error("vcftools_tajimad运行失败")
        self.option("stat_file", self.stat_path + ".Tajima.D")

    def run_vcftools_weir_fst_pop(self):
        """
        weir-fst-pop方法
        """
        self.stat_path = self.output_dir + "/" + os.path.basename(self.option("group_list")).split(".")[0] + "-"\
                         + os.path.basename(self.option("compare_group")).split(".")[0]
        cmd = "{} --vcf {} --weir-fst-pop {}".format(self.vcftools, self.option("vcf_file").prop["path"], self.option("group_list"))
        cmd += " --weir-fst-pop {} --fst-window-size {}".format(self.option("compare_group"), self.option("window_size"))
        cmd += " --fst-window-step {} --out {}".format(self.option("window_step"), self.stat_path)
        command = self.add_command("vcftools_weir_fst_pop", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("vcftools_weir_fst_pop运行成功")
        else:
            self.set_error("vcftools_weir_fst_pop运行失败")
        self.option("stat_file", self.stat_path + ".windowed.weir.fst")

    def run(self):
        super(VcftoolsStatTool, self).run()
        if self.option("analysis_method") == "window_pi":
            self.run_vcftools_window_pi()
        elif self.option("analysis_method") == "weir_fst_pop":
            self.run_vcftools_weir_fst_pop()
        else:
            self.run_vcftools_tajima_d()
        self.end()

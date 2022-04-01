# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# modified 2018.09.18

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os


class SweepRegionStatAgent(Agent):
    """
    工具： result.table.pl
    通过result.table.pl对结果进行统计
    """
    def __init__(self, parent):
        super(SweepRegionStatAgent, self).__init__(parent)
        options = [
            {"name": "vcf_file", "type": "infile", "format": "dna_gmap.vcf", "required": True},  # pop.recode.vcf或者比较分析的结果
            {"name": "pop_summary", "type": "infile", "format": "dna_gmap.vcf", "required": True},  # pop.summary
            {"name": "region_file", "type": "infile", "format": "dna_evolution.test", "required": True},  # pop1_pop2.1.pi_tajimaD_fst.select.region
            {"name": "window_weir_fst", "type": "infile", "format": "dna_evolution.window_weir_fst"},  # 1-2.windowed.weir.fst
            {"name": "tajima", "type": "infile", "format": "dna_evolution.tajima_d", "required": True},  # tajima_d
            {"name": "window_pi", "type": "infile", "format": "dna_evolution.window_pi", "required": True},  # window_pi
            {"name": "compare_tajima", "type": "infile", "format": "dna_evolution.tajima_d"},  # tajima_d
            {"name": "compare_window_pi", "type": "infile", "format": "dna_evolution.window_pi"},  # window_pi
            {"name": "compare_region_file", "type": "infile", "format": "dna_evolution.test"},  # pop1_pop2.2.pi_tajimaD_fst.select.region
            {"name": "pop_type", "type": "string", "default": "diff_pop"},  # 群体类型,diff_pop/pop
        ]
        self.add_option(options)

    def check_options(self):
        if self.option("pop_type") not in ["pop", "diff_pop"]:
            raise OptionError("pop_type:{}只能是diff_pop/pop，请检查")
        if self.option("compare_region_file").is_set and self.option("pop_type") != "diff_pop":
            raise OptionError("设置了compare_region_file文件的时候pop_type必须为diff_pop")
        if self.option("pop_type") == "diff_pop":
            if not self.option("compare_window_pi").is_set:
                raise OptionError("pop_type为diff_pop的时候必须设置compare_window_pi文件")
            if not self.option("compare_tajima").is_set:
                raise OptionError("pop_type为diff_pop的时候必须设置compare_tajima文件")
            if not self.option("window_weir_fst").is_set:
                raise OptionError("pop_type为diff_pop的时候必须设置window_weir_fst文件")

    def set_resource(self):
        self._cpu = 3
        self._memory = "15G"

    def end(self):
        super(SweepRegionStatAgent, self).end()


class SweepRegionStatTool(Tool):
    def __init__(self, config):
        super(SweepRegionStatTool, self).__init__(config)
        self.perl = "program/perl/perls/perl-5.24.0/bin/perl"
        self.perl_path = self.config.SOFTWARE_DIR + "/program/perl/perls/perl-5.24.0/bin/perl"
        self.parafly = "/program/parafly-r2013-01-21/src/ParaFly"
        self.result_table_path = self.config.PACKAGE_DIR + "/dna_evolution/result.table.pl"
        self.single_result_table_path = self.config.PACKAGE_DIR + "/dna_evolution/result.single_pop.table.pl"

    def run_result_table_single_pop(self):
        """
        单个群体的选择性消除结果统计：result.single_pop.table.pl
        """
        anno_stat = os.path.join(self.output_dir, os.path.basename(self.option("region_file").prop["path"]) + ".stat")
        cmd = "{} {} -vcf {}".format(self.perl, self.single_result_table_path, self.option("vcf_file").prop["path"])
        cmd += " -ann {} -region {}".format(self.option("pop_summary").prop["path"], self.option("region_file").prop["path"])
        cmd += " -pi {} -tajimad {}".format(self.option("window_pi").prop["path"], self.option("tajima").prop["path"])
        cmd += " -out {}".format(anno_stat)
        command = self.add_command("result_table_single_pop", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("result.single_pop.table.pl运行成功")
        else:
            self.set_error("result.single_pop.table.pl运行失败")

    def run_result_table_single(self):
        """
        result.table.pl
        """
        anno_stat = os.path.join(self.output_dir, os.path.basename(self.option("region_file").prop["path"]) + ".stat")
        cmd = "{} {} -vcf {}".format(self.perl, self.result_table_path, self.option("vcf_file").prop["path"])
        cmd += " -region {} -fst {}".format(self.option("region_file").prop["path"], self.option("window_weir_fst").prop["path"])
        cmd += " -pi1 {} -pi2 {}".format(self.option("window_pi").prop["path"], self.option("compare_window_pi").prop["path"])
        cmd += " -tajimad1 {} -tajimad2 {}".format(self.option("tajima").prop["path"], self.option("compare_tajima").prop["path"])
        cmd += " -ann {} -out {}".format(self.option("pop_summary").prop["path"], anno_stat)
        command = self.add_command("result_table_single", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("result_table_single运行成功")
        else:
            self.set_error("result_table_single运行失败")

    def run_result_table_double(self):
        """
        result.table.pl
        """
        anno_stat1 = os.path.join(self.output_dir, os.path.basename(self.option("region_file").prop["path"]) + ".stat")
        anno_stat2 = os.path.join(self.output_dir, os.path.basename(self.option("compare_region_file").prop["path"]) + ".stat")
        cmd1 = "{} {} -vcf {}".format(self.perl_path, self.result_table_path, self.option("vcf_file").prop["path"])
        cmd1 += " -region {} -fst {}".format(self.option("region_file").prop["path"], self.option("window_weir_fst").prop["path"])
        cmd1 += " -pi1 {} -pi2 {}".format(self.option("window_pi").prop["path"], self.option("compare_window_pi").prop["path"])
        cmd1 += " -tajimad1 {} -tajimad2 {}".format(self.option("tajima").prop["path"], self.option("compare_tajima").prop["path"])
        cmd1 += " -ann {} -out {}".format(self.option("pop_summary").prop["path"], anno_stat1)
        cmd2 = "{} {} -vcf {}".format(self.perl_path, self.result_table_path, self.option("vcf_file").prop["path"])
        cmd2 += " -region {} -fst {}".format(self.option("compare_region_file").prop["path"], self.option("window_weir_fst").prop["path"])
        cmd2 += " -pi1 {} -pi2 {}".format(self.option("window_pi").prop["path"], self.option("compare_window_pi").prop["path"])
        cmd2 += " -tajimad1 {} -tajimad2 {}".format(self.option("tajima").prop["path"], self.option("compare_tajima").prop["path"])
        cmd2 += " -ann {} -out {}".format(self.option("pop_summary").prop["path"], anno_stat2)
        cmd_file = os.path.join(self.work_dir, "result_table_cmd.list")
        wrong_cmd = os.path.join(self.work_dir, "result_table_cmd.txt")
        with open(cmd_file, "w") as f:
            f.write(cmd1 + "\n")
            f.write(cmd2 + "\n")
        cmd_more = "{} -c {} -CPU {} -failed_cmds {}".format(self.parafly, cmd_file, 2, wrong_cmd)
        command = self.add_command("result_table_double", cmd_more).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("{}运行完成".format("result_table_double"))
        else:
            self.set_error("{}运行失败".format("result_table_double"))

    def run(self):
        super(SweepRegionStatTool, self).run()
        if self.option("compare_region_file").is_set:
            self.run_result_table_double()
        elif self.option("pop_type") == "pop":
            self.run_result_table_single_pop()
        else:
            self.run_result_table_single()
        self.end()

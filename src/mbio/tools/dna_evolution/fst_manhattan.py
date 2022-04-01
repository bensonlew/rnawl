# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# modified 2018.08.22

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os


class FstManhattanAgent(Agent):
    """
    工具： fst-pi.R、 manhattan.R
    对vcf window-pi、TajimaD、weir-fst-pop出来的结果进行处理
    """
    def __init__(self, parent):
        super(FstManhattanAgent, self).__init__(parent)
        options = [
            {"name": "tajima", "type": "infile", "format": "dna_evolution.tajima_d"},  # tajima_d
            {"name": "window_pi", "type": "infile", "format": "dna_evolution.window_pi", "required": True},  # window_pi
            {"name": "compare_tajima", "type": "infile", "format": "dna_evolution.tajima_d"},  # tajima_d
            {"name": "compare_window_pi", "type": "infile", "format": "dna_evolution.window_pi", "required": True},  # window_pi
            {"name": "window_weir_fst", "type": "infile", "format": "dna_evolution.window_weir_fst", "required": True},  # window_weir_fst
            {"name": "p_value", "type": "float", "default": 0.05},  # p值筛选
            {"name": "analysis_type", "type": "string", "default": "all"},  # ["manhattan", "fst_pi", "all"]
            {"name": "select_file", "type": "outfile", "format": "dna_evolution.region_snp"},  # pop1_pop2.1.pi_tajimaD_fst.select
            {"name": "compare_select_file", "type": "outfile", "format": "dna_evolution.region_snp"},  # pop1_pop2.2.pi_tajimaD_fst.select
        ]
        self.add_option(options)

    def check_options(self):
        if self.option("p_value") <= 0 or self.option("p_value") >= 1:
            raise OptionError("p_value:%s的范围是0-1,请检查" % self.option("p_value"))
        if self.option("analysis_type") not in ["manhattan", "fst_pi", "all"]:
            raise OptionError("analysis_type:%s只能是manhattan/fst_pi/all,请检查" % self.option("analysis_type"))
        if self.option("analysis_type") != "fst_pi":
            if not self.option("tajima").is_set or not self.option("compare_tajima").is_set:
                raise OptionError("亲设置tajima文件或者compare_tajima文件")

    def set_resource(self):
        self._cpu = 2
        self._memory = "15G"

    def end(self):
        super(FstManhattanAgent, self).end()


class FstManhattanTool(Tool):
    def __init__(self, config):
        super(FstManhattanTool, self).__init__(config)
        self.R_path = "program/R-3.3.3/bin/Rscript"
        self.fst_pi_path = self.config.PACKAGE_DIR + "/dna_evolution/fst-pi.R"
        self.manhattan_path = self.config.PACKAGE_DIR + "/dna_evolution/manhattan.R"

    def run_fst_pi(self):
        """
        fst-pi.R
        """
        stat_path = self.output_dir + "/" + os.path.basename(self.option("window_weir_fst").prop["path"]).split(".")[0]
        cmd = "{} {} --fst {}".format(self.R_path, self.fst_pi_path, self.option("window_weir_fst").prop["path"])
        cmd += " --pi1 {} --pi2 {}".format(self.option("window_pi").prop["path"], self.option("compare_window_pi").prop["path"])
        cmd += " --thre {} --out {}".format(self.option("p_value"), stat_path)
        command = self.add_command("fst_pi", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("fst_pi运行成功")
        else:
            self.set_error("fst_pi运行失败")
        self.option("select_file", stat_path + ".fst_pi.detail.select")

    def run_manhattan(self):
        """
        manhattan.R
        """
        stat_path = self.output_dir + "/" + os.path.basename(self.option("window_weir_fst").prop["path"]).split(".")[0]
        cmd = "{} {} --fst {}".format(self.R_path, self.manhattan_path, self.option("window_weir_fst").prop["path"])
        cmd += " --pi1 {} --pi2 {}".format(self.option("window_pi").prop["path"], self.option("compare_window_pi").prop["path"])
        cmd += " --tajima1 {} --tajima2 {}".format(self.option("tajima").prop["path"], self.option("compare_tajima").prop["path"])
        cmd += " --thre {} --out {}".format(self.option("p_value"), stat_path)
        command = self.add_command("manhattan", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("manhattan运行成功")
        else:
            self.set_error("manhattan运行失败")
        self.option("select_file", stat_path + ".1.pi_tajimaD_fst.select")
        self.option("compare_select_file", stat_path + ".2.pi_tajimaD_fst.select")

    def run(self):
        super(FstManhattanTool, self).run()
        if self.option("window_weir_fst").prop["null"]:
            self.logger.info("{}文件为空，不进行fst计算".format(self.option("window_weir_fst").prop["path"]))
            self.end()
        if self.option("analysis_type") == "manhattan":
            self.run_manhattan()
        elif self.option("analysis_type") == "fst_pi":
            self.run_fst_pi()
        else:
            self.run_fst_pi()
            self.run_manhattan()
        self.end()

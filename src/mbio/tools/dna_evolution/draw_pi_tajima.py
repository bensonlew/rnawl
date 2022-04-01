# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# modified 2018.08.22

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os


class DrawPiTajimaAgent(Agent):
    """
    工具：pi-tajima.R
    对进行vcftools window-pi、TajimaD分析的结果进行处理
    """
    def __init__(self, parent):
        super(DrawPiTajimaAgent, self).__init__(parent)
        options = [
            {"name": "tajima_file", "type": "infile", "format": "dna_evolution.tajima_d", "required": True},  # tajima_d
            {"name": "window_pi_file", "type": "infile", "format": "dna_evolution.window_pi", "required": True},  # window_pi
            {"name": "p_value", "type": "float", "default": 0.05},  # p值筛选
            {"name": "select_file", "type": "outfile", "format": "dna_evolution.region_snp"},  # 4.pi_tajimaD.detail.select
        ]
        self.add_option(options)

    def check_options(self):
        if self.option("p_value") <= 0 or self.option("p_value") >= 1:
            raise OptionError("p_value:%s的范围是0-1,请检查" % self.option("p_value"))

    def set_resource(self):
        self._cpu = 2
        self._memory = "5G"

    def end(self):
        super(DrawPiTajimaAgent, self).end()


class DrawPiTajimaTool(Tool):
    def __init__(self, config):
        super(DrawPiTajimaTool, self).__init__(config)
        self.R_path = "program/R-3.3.3/bin/Rscript"
        self.pi_tajima_path = self.config.PACKAGE_DIR + "/dna_evolution/pi-tajima.R"

    def run_pi_tajima(self):
        """
        window-pi方法
        """
        self.stat_path = self.output_dir + "/" + os.path.basename(self.option("tajima_file").prop["path"]).split(".")[0]
        cmd = "{} {} --tajima {}".format(self.R_path, self.pi_tajima_path, self.option("tajima_file").prop["path"])
        cmd += " --pi {} --thre {} --out {}".format(self.option("window_pi_file").prop["path"], self.option("p_value"), self.stat_path)
        command = self.add_command("pi_tajima", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("pi_tajima运行成功")
        else:
            self.set_error("pi_tajima运行失败")
        self.option("select_file", self.stat_path + ".pi_tajimaD.detail.select")

    def run(self):
        super(DrawPiTajimaTool, self).run()
        if self.option("tajima_file").prop["null"]:
            self.logger.info("{}文件为空，不进行pi_tajimad计算".format(self.option("tajima_file").prop["path"]))
        elif self.option("window_pi_file").prop["null"]:
            self.logger.info("{}文件为空，不进行pi_tajimad计算".format(self.option("window_pi_file").prop["path"]))
        else:
            self.run_pi_tajima()
        self.end()

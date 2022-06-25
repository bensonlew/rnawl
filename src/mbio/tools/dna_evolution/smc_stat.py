# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# modified 2018.08.29

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import glob


class SmcStatAgent(Agent):
    """
    工具：drawPSMC.R对psmc结果进行整
    """
    def __init__(self, parent):
        super(SmcStatAgent, self).__init__(parent)
        options = [
            {"name": "smc_dir", "type": "infile", "format": "dna_evolution.sweep_dir", "required": True},
        ]
        self.add_option(options)

    def check_options(self):
        pass

    def set_resource(self):
        self._cpu = 2
        self._memory = "15G"

    def end(self):
        super(SmcStatAgent, self).end()


class SmcStatTool(Tool):
    def __init__(self, config):
        super(SmcStatTool, self).__init__(config)
        self.R_path = "program/R-3.3.3/bin/Rscript"
        self.perl = "miniconda2/bin/perl"
        self.draw_psmc = self.config.PACKAGE_DIR + "/dna_evolution/drawPSMC.R"
        self.merge_psmc = self.config.PACKAGE_DIR + "/dna_evolution/mergePSMCresult.pl"

    def run_all_smc_stat(self):
        """
        drawPSMC.R
        """
        cmd = "{} {} -in {} -out {}".format(self.perl, self.merge_psmc, self.option("smc_dir").prop["path"], os.path.join(self.output_dir, "all_smc_stat.xls"))
        command = self.add_command("all_smc_stat", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("all_smc_stat运行成功")
        else:
            self.set_error("all_smc_stat运行失败")

    # def set_output(self):
    #     files = glob.glob(r"*.smc.result")
    #     for f in files:
    #         f1 = os.path.join(self.work_dir, f)
    #         f2 = os.path.join(self.output_dir, f)
    #         if os.path.exists(f2):
    #             os.remove(f2)
    #         os.link(f1, f2)

    def run(self):
        super(SmcStatTool, self).run()
        self.run_all_smc_stat()
        self.end()

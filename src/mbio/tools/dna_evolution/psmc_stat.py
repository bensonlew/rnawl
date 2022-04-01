# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# modified 2018.08.29

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import glob


class PsmcStatAgent(Agent):
    """
    工具：drawPSMC.R对psmc结果进行整
    """
    def __init__(self, parent):
        super(PsmcStatAgent, self).__init__(parent)
        options = [
            {"name": "psmc_list", "type": "infile", "format": "dna_evolution.group_table", "required": True},  # psmc.list，组\tpsmc文件
        ]
        self.add_option(options)

    def check_options(self):
        pass

    def set_resource(self):
        self._cpu = 2
        self._memory = "15G"

    def end(self):
        super(PsmcStatAgent, self).end()


class PsmcStatTool(Tool):
    def __init__(self, config):
        super(PsmcStatTool, self).__init__(config)
        self.R_path = "program/R-3.3.3/bin/Rscript"
        self.perl = "program/perl/perls/perl-5.24.0/bin/perl"
        self.draw_psmc = self.config.PACKAGE_DIR + "/dna_evolution/drawPSMC.R"
        self.merge_psmc = self.config.PACKAGE_DIR + "/dna_evolution/mergePSMCresult.pl"

    def get_group_psmc(self):
        """
        将每个分组对应样本的psmc结果cat到一起，得到分组的psmc结果
        """
        group_dict = {}
        with open(self.option("psmc_list").prop["path"], "r") as f:
            lines = f.readlines()
            for line in lines[1:]:
                item = line.strip().split(" ")
                if item[0] not in group_dict.keys():
                    group_dict[item[0]] = []
                group_dict[item[0]].append(item[1])
        self.group_psmc_list = os.path.join(self.work_dir, "group_psmc.list")
        with open(self.group_psmc_list, "w") as w:
            w.write("popid\tfile\n")
            for group in group_dict.keys():
                psmc_group = os.path.join(self.work_dir, group + ".psmc")
                os.system("cat {} > {}".format(" ".join(group_dict[group]), psmc_group))
                w.write(group + "\t" + psmc_group + "\n")

    def run_draw_psmc(self):
        """
        drawPSMC.R
        """
        cmd = "{} {} --infile {} --outdir {}".format(self.R_path, self.draw_psmc, self.group_psmc_list, self.work_dir)
        command = self.add_command("draw_psmc", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("draw_psmc运行成功")
        else:
            self.set_error("draw_psmc运行失败")

    def run_all_psmc_stat(self):
        """
        drawPSMC.R
        """
        cmd = "{} {} -in {} -out {}".format(self.perl, self.merge_psmc, self.output_dir, os.path.join(self.output_dir, "all_psmc_stat.xls"))
        command = self.add_command("all_psmc_stat", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("all_psmc_stat运行成功")
        else:
            self.set_error("all_psmc_stat运行失败")

    def set_output(self):
        files = glob.glob(r"*.psmc.result")
        for f in files:
            f1 = os.path.join(self.work_dir, f)
            f2 = os.path.join(self.output_dir, f)
            if os.path.exists(f2):
                os.remove(f2)
            os.link(f1, f2)

    def run(self):
        super(PsmcStatTool, self).run()
        self.get_group_psmc()
        self.run_draw_psmc()
        self.set_output()
        self.run_all_psmc_stat()
        self.end()

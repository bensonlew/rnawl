# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# modified 2018.08.29

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os


class PsmcAgent(Agent):
    """
    工具：psmc
    """
    def __init__(self, parent):
        super(PsmcAgent, self).__init__(parent)
        options = [
            {"name": "psmcfa_list", "type": "infile", "format": "dna_evolution.group_table", "required": True},  # 每个组对应的psmcfa文件
            {"name": "iteration_num", "type": "int", "default": 25},  # -N,最大迭代次数
            {"name": "coalescent_time", "type": "float", "default": 25},  # -t,最大2NO合并时间
            {"name": "theta_ratio", "type": "float", "default": 5},  # -r,初始的theta/rho比率
            {"name": "parameters_pattern", "type": "string", "default": "4+25*2+4+6"},  # -p,参数模式
            {"name": "psmc_list", "type": "outfile", "format": "dna_evolution.group_table"},  # psmc.list
        ]
        self.add_option(options)

    def check_options(self):
        pass

    def set_resource(self):
        self.logger.info(self.option("psmcfa_list").prop["group_info"])
        group_num = len(self.option("psmcfa_list").prop["group_info"].keys())
        self._cpu = group_num + 1
        self._memory = "15G"

    def end(self):
        super(PsmcAgent, self).end()


class PsmcTool(Tool):
    def __init__(self, config):
        super(PsmcTool, self).__init__(config)
        self.psmc = self.config.SOFTWARE_DIR + "/bioinfo/dna_evolution/psmc"
        self.parafly = "/program/parafly-r2013-01-21/src/ParaFly"

    def run_psmc(self):
        """
        psmc
        """
        self.psmcfa_info = {}
        cmd_list = []
        with open(self.option("psmcfa_list").prop["path"], "r") as f, open(self.output_dir + "/psmc.list", "w") as w:
            w.write("popid file\n")
            for line in f:
                item = line.strip().split("\t")
                self.psmcfa_info[item[0]] = item[1]
                w.write(item[0] + " " + os.path.join(self.output_dir, item[0] + ".psmc") + "\n")
                # w.write("all" + " " + os.path.join(self.output_dir, item[0] + ".psmc") + "\n")
        for name in self.psmcfa_info.keys():
            psmc_file = os.path.join(self.output_dir, name + ".psmc")
            cmd = "{} -N{} -t{}".format(self.psmc, self.option("iteration_num"), self.option("coalescent_time"))
            cmd += " -r{} -p {} -o {}".format(self.option("theta_ratio"), self.option("parameters_pattern"), psmc_file)
            cmd += " {}".format(self.psmcfa_info[name])
            cmd_list.append(cmd)
        cmd_file = os.path.join(self.work_dir, "psmc_cmd.list")
        wrong_cmd = os.path.join(self.work_dir, "failed_psmc_cmd.txt")
        with open(cmd_file, "w") as f:
            for cmd in cmd_list:
                f.write(cmd + "\n")
        cmd_more = "{} -c {} -CPU {} -failed_cmds {}".format(self.parafly, cmd_file, len(cmd_list), wrong_cmd)
        command = self.add_command("psmc_more", cmd_more).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("{}运行完成".format("psmc_more"))
        else:
            self.set_error("{}运行失败".format("psmc_more"))
        self.option("psmc_list", self.output_dir + "/psmc.list")

    def run(self):
        super(PsmcTool, self).run()
        self.run_psmc()
        self.end()

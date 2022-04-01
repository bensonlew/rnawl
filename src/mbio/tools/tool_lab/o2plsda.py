# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
# version 1.0
# last_modify: 2021.05.08

import os, time
import re, shutil
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from mbio.packages.bac_comp_genome.common_function import link_dir,link_file
import json
import subprocess

class O2plsdaAgent(Agent):
    """
    代谢O2plsda分析
    """
    def __init__(self, parent):
        super(O2plsdaAgent, self).__init__(parent)
        options = [
            {"name": "x_list", "type": "string"},
            {"name": "y_list", "type": "string"},
            {"name": "x_data", "type": "infile", "format": "sequence.profile_table"},
            {"name": "y_data", "type": "infile", "format": "sequence.profile_table"},
            {"name": "group", "type": "infile", "format": "sequence.profile_table" },
            {"name": "oxoy", "type": "string", 'default':"Microbe;Metabolome"},
            {"name": "scale", "type": "string", 'default': "UV"},
            {"name": "x_method", "type": "string", 'default': "log1"},
            {"name": "y_method", "type": "string", 'default': "log1"},
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("x_data").is_set:
            raise OptionError("必须输入x_data文件")

    def set_resource(self):
        self._cpu = 2
        self._memory = '10G'

    def end(self):
        super(O2plsdaAgent, self).end()

class O2plsdaTool(Tool):
    def __init__(self, config):
        super(O2plsdaTool, self).__init__(config)
        self.path = self.config.SOFTWARE_DIR + "/program/R-3.3.1/bin"
        self.lib = self.config.SOFTWARE_DIR + "/program/R-3.3.3/lib64/R/lib"
        self.set_environ(PATH=self.path, LD_LIBRARY_PATH=self.lib)
        self.python = "/program/Python/bin/python"
        self.R = "/program/R-3.3.1/bin/Rscript"
        self.script = self.config.PACKAGE_DIR + '/tool_lab/'

    def run_O2PLS(self):
        """
        运行O2PLS
        :return:
        """
        cmd = "{} {}run_O2PLS_TwoOmics.py -omics_x {} -omics_y {} -group {} -oxoy {} -scale {} -log_x {} -log_y {} -out {} -path {}".format(self.python, self.script, self.option("x_data").prop['path'], self.option("y_data").prop['path'],self.option("group").prop['path'],self.option("oxoy"), self.option("scale"), self.option("x_method"),self.option("y_method"),self.work_dir+"/result",self.script + "scaleNorm.r")
        self.logger.info(cmd)
        command = self.add_command("run_o2pls", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info('生成result.O2PLS.cmd.r成功！')
        else:
            self.set_error("run_O2PLS运行出错!")

    def run_r(self):
        cmd = "/program/R-3.3.1/bin/Rscript {}".format(self.work_dir+"/result.O2PLS.cmd.r")
        self.logger.info(cmd)
        command = self.add_command("run_r", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info('生成run_r成功！')
        else:
            self.set_error("run_r运行出错!")

    def run_rename(self):
        with open(self.option('x_list'),"r") as f,open(self.option('y_list'),"r") as s,open(self.work_dir + "/result.O2PLS_Loadings.xls","r") as k,open(self.output_dir + "/O2PLS_Loadings.xls","w") as r:
            dict ={}
            lines = f.readlines()
            for line in lines:
                lin =line.strip().split("\t")
                dict[lin[0]] = lin[1]
            lines2 = s.readlines()
            for line in lines2:
                lin = line.strip().split("\t")
                dict[lin[0]] = lin[1]
            lines3 = k.readlines()
            r.write(lines3[0])
            for line in lines3[1:]:
                lin = line.strip().split("\t")
                dd = dict[lin[0]]
                lin[0] = dd
                r.write("\t".join(lin)+"\n")

    def set_output(self):
        """
        设置结果文件目录
        """
        link_file(self.work_dir+ "/result.fit_summary.xls", self.output_dir+"/fit_summary.xls")
        link_file(self.work_dir + "/result.O2PLS_Scores.xls", self.output_dir + "/O2PLS_Scores.xls")
        self.logger.info("生成结果文件完成")

    def run(self):
        """
        运行
        """
        super(O2plsdaTool, self).run()
        self.run_O2PLS()
        self.run_r()
        self.run_rename()
        self.set_output()
        self.end()
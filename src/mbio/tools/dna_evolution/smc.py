# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# modified 2018.09.08

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import shutil


class SmcAgent(Agent):
    """
    工具：smc++ estimate/plot/结果重写
    得到一个分群的smc结果，用于多个分群统一
    """
    def __init__(self, parent):
        super(SmcAgent, self).__init__(parent)
        options = [
            {"name": "smc_dir", "type": "infile", "format": "dna_evolution.sweep_dir", "required": True},  # 一个分群的所有的smc结果，不能是空文件夹
        ]
        self.add_option(options)

    def check_options(self):
        pass

    def set_resource(self):
        self._cpu = 2
        self._memory = "50G"

    def end(self):
        super(SmcAgent, self).end()


class SmcTool(Tool):
    def __init__(self, config):
        super(SmcTool, self).__init__(config)
        self.set_environ(PATH=self.config.SOFTWARE_DIR + '/gcc/7.2.0/bin')
        self.set_environ(LD_LIBRARY_PATH=self.config.SOFTWARE_DIR + '/gcc/7.2.0/lib64')
        self.set_environ(PATH=self.config.SOFTWARE_DIR + '/program/Python35/bin')
        self.set_environ(LD_LIBRARY_PATH=self.config.SOFTWARE_DIR + '/program/Python35/lib')
        self.set_environ(LD_LIBRARY_PATH=self.config.SOFTWARE_DIR + '/library/gcc_7.2/gmp-6.1.0/lib')
        self.set_environ(LD_LIBRARY_PATH=self.config.SOFTWARE_DIR + '/library/gcc_7.2/gsl23/lib')
        self.set_environ(LD_LIBRARY_PATH=self.config.SOFTWARE_DIR + '/library/gcc_7.2/mpfr-3.1.4/lib')
        self.smc_path = "program/Python35/bin/smc++"

    def run_smc_estimate(self):
        """
        smc++ estimate 方法
        smc++ estimate -o ./ 1.25e-8 *.smc.gz
        """
        # self.logger.info(os.environ)
        smc_dir = []
        for f in os.listdir(self.option("smc_dir").prop["path"]):
            if f.endswith(".smc.gz"):
                smc_dir.append(os.path.join(self.option("smc_dir").prop["path"], f))
        self.logger.info(smc_dir)
        cmd = "{} estimate -o {} 1.25e-8 {}".format(self.smc_path, self.output_dir, ' '.join(smc_dir))
        command = self.add_command("smc_estimate", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("smc_estimate运行成功")
        else:
            self.set_error("smc_estimate运行失败")

    def run_smc_plot(self):
        """
        smc++ plot 方法
        smc++ plot --csv ./test model.final.json
        """
        pop_name = os.path.basename(self.option("smc_dir").prop["path"])
        cmd = "{} plot --csv {} ".format(self.smc_path, os.path.join(self.output_dir, pop_name))
        cmd += os.path.join(self.output_dir, "model.final.json")
        command = self.add_command("smc_plot", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("smc_plot运行成功")
        else:
            self.set_error("smc_plot运行失败")

    def run_smc_plot_deale(self):
        """
        smc++ plot结果处理，取x,y两列组成新文件
        """
        pop_name = os.path.basename(self.option("smc_dir").prop["path"])
        with open(os.path.join(self.output_dir, pop_name + ".csv"), "r") as f, open(os.path.join(self.output_dir, pop_name + ".smc.result"), "w") as w:
            lines = f.readlines()
            w.write("\"Generation x 1\"\t\"Ne(log 10)\"\n")
            for line in lines[1:]:
                item = line.strip().split(",")
                w.write(item[1] + "\t" + item[2] + "\n")

    def run(self):
        super(SmcTool, self).run()
        self.run_smc_estimate()
        self.run_smc_plot()
        self.run_smc_plot_deale()
        self.end()

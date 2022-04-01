# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
# version 1.0
# last_modify: 2021.04.29

import os, time
import re, shutil
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from mbio.packages.bac_comp_genome.common_function import link_dir
import json
import pandas as pd

class SerotypefinderAgent(Agent):
    """
    单个基因组serotypefinder预测
    """
    def __init__(self, parent):
        super(SerotypefinderAgent, self).__init__(parent)
        options = [
            {"name": "fasta", "type": "infile", "format": "sequence.fasta"},# 序列文件
            {"name": "sample_name", "type": "string"},# 样本名
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("fasta").is_set:
            raise OptionError("必须输入序列文件")

    def set_resource(self):
        self._cpu = 2
        self._memory = '20G'

    def end(self):
        super(SerotypefinderAgent, self).end()

class SerotypefinderTool(Tool):
    def __init__(self, config):
        super(SerotypefinderTool, self).__init__(config)
        self.path = self.config.SOFTWARE_DIR + "/program/Python35/bin"
        self.lib = self.config.SOFTWARE_DIR + "/program/Python35/lib"
        self.set_environ(PATH=self.path, LD_LIBRARY_PATH=self.lib)
        self.python = "/program/Python/bin/python"
        self.python2 = "/program/Python35/bin/python"
        self.python_script = self.config.SOFTWARE_DIR + "/bioinfo/Genomic/Sofware/serotypefinder/"
        self.blastn = self.config.SOFTWARE_DIR + "/bioinfo/align/ncbi-blast-2.3.0+/bin/blastn"
        self.kma = self.config.SOFTWARE_DIR + "/bioinfo/Genomic/Sofware/kma/"
        self.db = self.config.SOFTWARE_DIR + "/bioinfo/Genomic/Sofware/serotypefinder/serotypefinder_db/"


    def run_serotypefinder(self):
        """
        运行mlst
        :return:
        """
        if os.path.exists(self.work_dir+"/result"):
           shutil.rmtree(self.work_dir+"/result")
        os.mkdir(self.work_dir+"/result")
        cmd = "{} {}serotypefinder.py -i {} -o {} -mp {} -p {} -q -x".format(self.python2, self.python_script, self.option("fasta").prop['path'], self.work_dir+"/result", self.blastn, self.db)
        self.logger.info(cmd)
        command = self.add_command("run_serotypefinder", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_serotypefinder运行完成！")
        else:
            self.set_error("run_serotypefinder运行出错!")

    def run_stat(self):
        """
        运行处理文件，生成结果文件
        :return:
        """
        if os.path.exists(self.work_dir + "/result/data.json"):
            cmd = "{} {}serotypefinder_stat.py --d {} --s {} --o {}".format(self.python, self.python_script,
                                                                            self.work_dir + "/result/data.json",
                                                                            self.option("sample_name"),
                                                                            self.output_dir + "/" + self.option(
                                                                                "sample_name"))
            self.logger.info(cmd)
            command = self.add_command("run_stat", cmd).run()
            self.wait(command)
            if command.return_code == 0:
                self.logger.info("run_stat运行完成！")
            else:
                self.set_error("run_stat运行出错!")


    def set_output(self):
        """
        设置结果文件目录
        """
        os.link(self.work_dir + "/result/Hit_in_genome_seq.fsa",
                self.output_dir + "/" + self.option("sample_name") + ".HitMLST.fasta")
        self.logger.info("生成结果文件完成")


    def run(self):
        """
        运行
        """
        super(SerotypefinderTool, self).run()
        self.run_serotypefinder()
        self.run_stat()
        self.set_output()
        self.end()
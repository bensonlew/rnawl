# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
# version 1.0
# last_modify: 2021.04.08

import os, time
import re, shutil
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from mbio.packages.bac_comp_genome.common_function import link_dir
import json

class CrisprcasfinderAgent(Agent):
    """
    单个基因组Crisprcasfinder预测
    """
    def __init__(self, parent):
        super(CrisprcasfinderAgent, self).__init__(parent)
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
        self._memory = '10G'

    def end(self):
        super(CrisprcasfinderAgent, self).end()

class CrisprcasfinderTool(Tool):
    def __init__(self, config):
        super(CrisprcasfinderTool, self).__init__(config)
        self.path = self.config.SOFTWARE_DIR + "/bioinfo/Genomic/Sofware/CRISPRCasFinder/bin:" +self.config.SOFTWARE_DIR + "/program/perl/perls/perl-5.24.0/bin"
        self.lib = self.config.SOFTWARE_DIR + "/program/perl/perls/perl-5.24.0/lib"
        self.set_environ(PATH=self.path, LD_LIBRARY_PATH=self.lib)
        self.db = self.config.SOFTWARE_DIR + "/bioinfo/Genomic/Sofware/CRISPRCasFinder/src/sel392v2.so"
        self.db2 = self.config.SOFTWARE_DIR + "/bioinfo/Genomic/Sofware/CRISPRCasFinder/CasFinder-2.0.3"
        self.perl = "/program/perl/perls/perl-5.24.0/bin/perl"
        self.python    = "/program/Python/bin/python"
        self.perl_script = self.config.SOFTWARE_DIR + "/bioinfo/Genomic/Sofware/CRISPRCasFinder/"
        self.script = self.config.PACKAGE_DIR + '/tool_lab/'

    def run_crisprcasfinder(self):
        """
        运行mlst
        :return:
        """
        if os.path.exists(self.work_dir+"/result"):
           shutil.rmtree(self.work_dir+"/result")
        cmd = "{} {}CRISPRCasFinder.pl -in {} -cf {} -def G -keep -so {} -out {}".format(self.perl, self.perl_script, self.option("fasta").prop['path'], self.db2, self.db, self.work_dir+"/result")
        self.logger.info(cmd)
        command = self.add_command("run_crisprcasfinder", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_crisprcasfinder运行完成！")
        else:
            self.set_error("run_crisprcasfinder运行出错!")

    def run_stat(self):
        """
        运行run_stat,主要统计crisprcasfinder的信息
        :return:
        """
        if os.path.exists(self.work_dir+"/result/result.json"):
            cmd = "{} {}crisprcasfinder_stat.py --d {} --s {} --o {}".format(self.python, self.script,
                                                                  self.work_dir + "/result/result.json",
                                                                  self.option("sample_name"),
                                                                  self.output_dir + "/" + self.option("sample_name"))
            self.logger.info(cmd)
            command = self.add_command("run_stat", cmd).run()
            self.wait(command)
            if command.return_code == 0:
                self.logger.info("run_stat运行完成！")
            else:
                self.set_error("run_stat运行出错!")
        else:
            self.end()


    def set_output(self):
        """
        设置结果文件目录
        """
        self.logger.info("生成结果文件完成")

    def run(self):
        """
        运行
        """
        super(CrisprcasfinderTool, self).run()
        self.run_crisprcasfinder()
        self.run_stat()
        self.set_output()
        self.end()
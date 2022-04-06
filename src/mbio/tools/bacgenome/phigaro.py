# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
# version 1.0
# last_modify: 2021.04.19

import os, time
import re, shutil
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from mbio.packages.bac_comp_genome.common_function import link_dir
import json

class PhigaroAgent(Agent):
    """
    单个基因组phigaro预测前噬菌体
    """
    def __init__(self, parent):
        super(PhigaroAgent, self).__init__(parent)
        options = [
            {"name": "fasta", "type": "infile", "format": "sequence.fasta"},# 序列文件
            {"name": "sample_name", "type": "string"},# 样本名
            {"name": "anno", "type": "infile", "format": "sequence.profile_table"},  # 注释总览表
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("fasta").is_set:
            raise OptionError("必须输入序列文件")

    def set_resource(self):
        self._cpu = 2
        self._memory = '10G'

    def end(self):
        super(PhigaroAgent, self).end()

class PhigaroTool(Tool):
    def __init__(self, config):
        super(PhigaroTool, self).__init__(config)
        self.path = self.config.SOFTWARE_DIR + "/program/Python35/bin"
        self.lib = self.config.SOFTWARE_DIR + "/program/Python35/lib"
        self.set_environ(PATH=self.path, LD_LIBRARY_PATH=self.lib)
        self.db = self.config.SOFTWARE_DIR + "/bioinfo/Genomic/Sofware/phigaro/config.yml"
        self.python = "/miniconda2/bin/python"
        self.python2 = "/program/Python35/bin/phigaro"
        self.script = self.config.PACKAGE_DIR + '/bacgenome/'

    def run_phigaro(self):
        """
        运行phigaro
        :return:
        """
        if os.path.exists(self.work_dir+"/"+self.option("sample_name")+".fasta"):
            os.remove(self.work_dir+"/"+self.option("sample_name")+".fasta")
        os.link(self.option("fasta").prop['path'],self.work_dir+"/"+self.option("sample_name")+".fasta")
        cmd = "{} -f {} -p -e tsv -o {} -c {} -d -t 4".format(self.python2, self.work_dir+"/"+self.option("sample_name")+".fasta", self.work_dir+"/result", self.db)
        self.logger.info(cmd)
        command = self.add_command("run_phigaro", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_phigaro运行完成！")
        else:
            self.set_error("run_phigaro运行出错!")
        with open(self.work_dir+"/result/"+self.option("sample_name")+".phigaro.tsv", "r") as f:
            lines = f.readlines()
            return len(lines)

    def run_stat(self):
        """
        运行run_stat,主要统计mlst的信息
        :return:
        """
        cmd = "{} {}phigaro_stat.py -i {} -d {} -f {} -p {}".format(self.python, self.script,
                                                                                self.work_dir+"/result/"+self.option("sample_name")+".phigaro.tsv",
                                                                                self.option("anno").prop['path'],self.option("fasta").prop['path'],
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
        self.logger.info("生成结果文件完成")

    def run(self):
        """
        运行
        """
        super(PhigaroTool, self).run()
        num = self.run_phigaro()
        if num<=1:
            self.end()
        else:
            self.run_stat()
            self.set_output()
            self.end()
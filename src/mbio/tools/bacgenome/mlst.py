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

class MlstAgent(Agent):
    """
    单个基因组MLST预测
    """
    def __init__(self, parent):
        super(MlstAgent, self).__init__(parent)
        options = [
            {"name": "fasta", "type": "infile", "format": "sequence.fasta"},# 序列文件
            {"name": "sample_name", "type": "string"},# 样本名
            {"name": "species", "type": "string", },# 物种名称
            {"name": "method", "type": "string", 'default':"blastn"},  #方法 kmer or blastn
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("fasta").is_set:
            raise OptionError("必须输入序列文件")

    def set_resource(self):
        self._cpu = 2
        self._memory = '10G'

    def end(self):
        super(MlstAgent, self).end()

class MlstTool(Tool):
    def __init__(self, config):
        super(MlstTool, self).__init__(config)
        self.path = self.config.SOFTWARE_DIR + "/program/Python35/bin" 
        self.lib = self.config.SOFTWARE_DIR + "/program/Python35/lib"
        self.set_environ(PATH=self.path, LD_LIBRARY_PATH=self.lib)
        self.db = self.config.SOFTWARE_DIR + "/bioinfo/Genomic/Sofware/mlst/mlst_database/"
        self.python = "/program/Python/bin/python"
        self.python2 = "/program/Python35/bin/python"
        self.python_script = self.config.SOFTWARE_DIR + "/bioinfo/Genomic/Sofware/mlst/mlst/"
        self.blastn = self.config.SOFTWARE_DIR + "/bioinfo/align/ncbi-blast-2.3.0+/bin/blastn"
        self.kma = self.config.SOFTWARE_DIR + "/bioinfo/Genomic/Sofware/kma"
        self.script = self.config.PACKAGE_DIR + '/bacgenome/'

    def run_mlst(self):
        """
        运行mlst
        :return:
        """
        if os.path.exists(self.work_dir+"/result"):
           shutil.rmtree(self.work_dir+"/result")
        os.mkdir(self.work_dir+"/result")
        if self.option('method') in ['blastn']:
            method = self.blastn
        elif self.option('method') in ['blastn']:
            method = self.kma
        cmd = "{} {}mlst.py -i {} -o {} -s {} -mp {} -x -p {}".format(self.python2, self.python_script, self.option("fasta").prop['path'], self.work_dir+"/result", self.option("species"), method, self.db)
        self.logger.info(cmd)
        command = self.add_command("run_mlst", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_mlst运行完成！")
        else:
            self.set_error("run_mlst运行出错!")

    def run_stat(self):
        """
        运行run_stat,主要统计mlst的信息
        :return:
        """
        cmd = "{} {}mlst_stat.py --d {} --s {} --o {}".format(self.python, self.script,
                                                                                self.work_dir+"/result/data.json",
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
        os.link(self.work_dir+"/result/Hit_in_genome_seq.fsa", self.output_dir+"/"+self.option("sample_name")+".HitMLST.fasta")
        self.logger.info("生成结果文件完成")

    def run(self):
        """
        运行
        """
        super(MlstTool, self).run()
        self.run_mlst()
        self.run_stat()
        self.set_output()
        self.end()
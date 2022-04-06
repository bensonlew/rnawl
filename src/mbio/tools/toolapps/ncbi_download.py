# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
# version 1.0
# last_modify: 2021.01.05

import os
import re
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import subprocess
import shutil


class NcbiDownloadAgent(Agent):
    """
    根据基因组NCBI编号下载数据
    """

    def __init__(self, parent):
        super(NcbiDownloadAgent, self).__init__(parent)
        options = [
            {"name": "sample_list", "type": "string"},  # 样品名称，sample1;sample2;...
            {"name": "genomes", "type": "outfile", "format": "sequence.fasta_dir"},
            {"name": "s16", "type": "outfile", "format": "sequence.fasta"},
            {"name": "table", "type": "outfile", "format": "sequence.profile_table"},
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("sample_list"):
            raise OptionError("必须设置参数sample_list!")

    def set_resource(self):
        self._cpu = 2
        self._memory = '10G'

    def end(self):
        super(NcbiDownloadAgent, self).end()


class NcbiDownloadTool(Tool):
    def __init__(self, config):
        super(NcbiDownloadTool, self).__init__(config)
        self.samplelist = self.option("sample_list")
        self.python = self.config.SOFTWARE_DIR + "/miniconda2/bin/python"
        self.script = self.config.PACKAGE_DIR + "/toolapps/"
        self.path = self.config.SOFTWARE_DIR + "/bioinfo/Genomic/Sofware/edirect:" +self.config.SOFTWARE_DIR + "/../.aspera/connect/bin:" +self.config.SOFTWARE_DIR + "/program/:"
        self.set_environ(PATH=self.path)

    def run_ncbidata(self):
        if os.path.exists(self.work_dir+"/ncbi"):
            shutil.rmtree(self.work_dir+"/ncbi")
        os.mkdir(self.work_dir+"/ncbi")
        cmd = "export PATH={}$PATH\n".format(self.path)
        cmd += "{} {}download_genome.py -ID_list \"{}\" -output {}".format(self.python,self.script, str(self.samplelist), self.work_dir+"/ncbi")
        try:
            self.logger.info(cmd)
            subprocess.check_output(cmd, shell=True)
            self.logger.info("download_genome运行完成")
        except subprocess.CalledProcessError:
            self.set_error("download_genome运行出错", code="33300201")

    def run_getfile(self):
        if os.path.getsize(self.work_dir+ "/file.txt") > 0:
            if os.path.exists(self.work_dir+"/genomes"):
                shutil.rmtree(self.work_dir+"/genomes")
            os.mkdir(self.work_dir+"/genomes")
            if os.path.exists(self.work_dir+"/16s"):
                shutil.rmtree(self.work_dir+"/16s")
            os.mkdir(self.work_dir + "/16s")
            des =[]
            with open(self.work_dir + "/file.txt","r") as f:
                lines = f.readlines()
                for line in lines:
                    lin = line.strip().split("\t")
                    os.link(lin[1], self.work_dir+"/genomes/"+lin[0]+".fasta")
                    if lin[3] != 0:
                        des.append(lin[4])
            if len(des) >0:
                os.system("cat {} >{}".format(" ".join(des), self.output_dir+"/all.database_16s.fasta"))

    def set_output(self):
        if os.path.exists(self.output_dir + "/file.txt"):
            os.remove(self.output_dir + "/file.txt")
        os.link(self.work_dir + "/file.txt", self.output_dir + "/file.txt")
        self.option("genomes", self.work_dir+"/genomes")
        self.option("s16", self.output_dir+"/all.database_16s.fasta")
        self.option("table", self.output_dir + "/file.txt")


    def run(self):
        super(NcbiDownloadTool, self).run()
        self.run_ncbidata()
        self.run_getfile()
        self.set_output()
        self.end()
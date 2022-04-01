# !usr/bin/python
# -*- coding: utf-8 -*-
# __author__ = 'XueQinwen'

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
#from biocluster.config import config
from collections import namedtuple, defaultdict
from Bio import SeqIO
import time
import re
import time
import random
import datetime
import subprocess
import unittest
import os
import glob
import sys
import shutil

class VcfDistributionV2Agent(Agent):
    def __init__(self, parent):
        super(VcfDistributionV2Agent, self).__init__(parent)
        options = [
            {"name":"vcf_file","type":"infile","format":"dna_gmap.vcf"},
            {"name":"wsize","type":"int","default":100000}
        ]
        self.add_option(options)
    
    def check_option(self):
        """
        参数设置
        """
        if not self.option("vcf_file"):
            raise OptionError("没有输入vcf文件")
    
    def set_resource(self):
        """
        设置所需资源
        """
        self._cpu = 2
        self._memory = "15G"

    def end(self):
        super(VcfDistributionV2Agent,self).end()


class VcfDistributionV2Tool(Tool):
    def __init__(self, config):
        super(VcfDistributionV2Tool, self).__init__(config)
        self._version = '2.0'
        # self.Rscript = os.path.join(self.config.SOFTWARE_DIR,"program/miniconda3/bin/activate")
        # self.Rscript = "program/miniconda3/envs/R_v4/bin/Rscript"
        self.Rscript = "bioinfo/ysource/miniconda3/envs/R_v4/bin/Rscript"
        self.Rmvp = self.config.PACKAGE_DIR + "/tool_lab/MVP.singleV2.R"
        self.shell_path = self.config.SOFTWARE_DIR+"/program/sh"
        self.set_environ(
            LD_LIBRARY_PATH=self.config.SOFTWARE_DIR + '/library/glibc-2.14/lib/')
        self.conda = self.config.SOFTWARE_DIR + "/program/miniconda3/etc/profile.d/conda.sh"
        
    def run(self):
        """
        运行
        """
        super(VcfDistributionV2Tool, self).run()
        self.run_Rscript()
        self.end()

    def run_Rscript(self):
        """
        Rscript MVP.singleV2.R --vcf vcffile --output output --wsize 100000
        """
        # now_time = datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f") + \
        #     "_" + str(random.randint(1, 10000))
        # script_path = self.work_dir
        # if not os.path.exists(script_path):
        #     os.mkdir(script_path)
        # file_path = script_path + "/distribution_{}.sh".format(now_time)
        # cmd = "source {}\n".format(self.conda)
        # cmd += "conda activate R_v4 \n"
        cmd = "{} {} --vcf {} ".format(self.Rscript,self.Rmvp, self.option("vcf_file").prop["path"])
        # cmd += "Rscript {} --plink {}/pop ".format(self.Rmvp, self.option("pop_dir").prop["path"])
        cmd += "--wsize {} ".format(self.option("wsize"))
        cmd += "--output {} \n".format(self.output_dir)
        # cmd += "conda deactivate\n"
    
        # self.logger.info(cmd)
        self.logger.info(cmd)
        command = self.add_command("rmvp", cmd).run()
        self.logger.info("开始画图")
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("画图完成")
        else:
            self.set_error("画图错误，请重新检查输入参数")
        # with open(file_path, 'w') as w:
        #     w.write('#!/bin/bash'+"\n")
        #     w.write(cmd)
        # code = os.system('chmod +x {}'.format(file_path))
        # if code == 0:
        #     self.logger.info("修改{}为可执行文件成功！".format(file_path))
        # else:
        #     self.set_error("修改{}为可执行文件失败！".format(file_path))
        # shell = "/bioinfo/ancestorage/script_temp/{}".format(
        #     os.path.basename(file_path))
        # cmd1 = "{} {}".format(self.shell_path, file_path)
        # self.logger.info(shell)
        # try:
        #     subprocess.check_output(cmd1, shell=True)
        # except subprocess.CalledProcessError:
        #     self.set_error("画图出现错误")
        # os.system('rm {}'.format(file_path))

    
    

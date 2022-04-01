# !usr/bin/python
# -*- coding: utf-8 -*-
# __author__ = 'wangwenjie'

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from collections import namedtuple, defaultdict
from Bio import SeqIO
import pandas as pd
import re
import unittest
import os
import glob
import sys
import shutil

class LmplotAgent(Agent):
    def __init__(self, parent):
        super(LmplotAgent, self).__init__(parent)
        options = [
            {"name": "input_file", "type": "infile", "format": "denovo_rna_v2.common"},
            {"name": "gene_name", "type": "string"}
        ]
        self.add_option(options)

    def check_option(self):
        """
        参数检查
        """
        if not self.option("input_file"):
            raise OptionError("没有找到input_file")
        if not self.option("gene_name"):
            raise OptionError("没有找到gene_name")
        return True

    def set_resource(self):
        """
        设置所需资源
        """
        self._cpu = 8
        self._memory = "5G"

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        super(LmplotAgent, self).end()

class LmplotTool(Tool):
    def __init__(self, config):
        super(LmplotTool, self).__init__(config)
        # self.Rscript = "bioinfo/ysource/miniconda2/envs/R4version/bin/Rscript"
        self.Rscript = "bioinfo/ysource/miniconda3/envs/R_v4/bin/Rscript"
        self.line_r = self.config.SOFTWARE_DIR+'/bioinfo/basic_molecular/line.R'

    def run(self):
        """
        运行
        """
        super(LmplotTool, self).run()
        self.run_lmplot()
        self.end()

    def run_lmplot(self):
        """
        计算标准曲线
        """
        gene_name = self.option("gene_name")
        cmd = "{} {} {} {}".format(self.Rscript, self.line_r, self.option("input_file").prop["path"], gene_name)
        # os.system('cd {} && mkdir lmplot_dir'.format(self.work_dir))
        self.logger.info(cmd)
        command = self.add_command("lmplot", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("计算标准曲线成功")
        else:
            self.logger.info("计算标准曲线失败")
        self.set_output()
        
    def set_output(self):
        gene_name = self.option("gene_name")
        os.link(os.path.join(self.work_dir, "{}_line.png".format(gene_name)),
                    os.path.join(self.output_dir, "{}_line.png".format(gene_name)))
        os.link(os.path.join(self.work_dir, "{}_lineOut.txt".format(gene_name)),
                    os.path.join(self.output_dir, "{}_lineOut.txt".format(gene_name)))

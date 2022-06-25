#!/usr/bin/env python
# -*- coding: utf-8 -*-
# __author__ = 'gaohao'

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import pandas as pd

class CheckmTraAgent(Agent):
    """
    用于评估bin和genome的完整性
    version 1.0
    author: gaohao
    last_modify: 2019.01.07
    """

    def __init__(self, parent):
        super(CheckmTraAgent, self).__init__(parent)
        options = [
            {"name": "input_genome", "type": "infile", "format": "sequence.fasta"},
            {"name": "out", "type": "outfile", "format": "meta.otu.otu_table, meta.otu.tax_summary_dir"},  # 生成结果文件
        ]
        self.add_option(options)

    def check_options(self):
        """
        检测参数是否正确
        """
        if not self.option('input_genome').is_set:
            raise OptionError("input_genome文件不存在！")

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 2
        self._memory = '10G'

    def end(self):
        super(CheckmTraAgent, self).end()


class CheckmTraTool(Tool):
    """
    version 1.0
    """
    def __init__(self, config):
        super(CheckmTraTool, self).__init__(config)
        self.path =self.config.SOFTWARE_DIR + "/bioinfo/metaGenomic/CheckM-master:" + self.config.SOFTWARE_DIR + "/bioinfo/metaGenomic/Prodigal-2.6.3:" + self.config.SOFTWARE_DIR + "/bioinfo/Genomic/Sofware/hmmer-3.1b2/bin:" + self.config.SOFTWARE_DIR + "/bioinfo/metaGenomic/pplacer-Linux-v1.1.alpha19"
        self.set_environ(PATH=self.path)
        self.checkm ="/miniconda2/bin/"
        self.perl = "/miniconda2/bin/perl"
        self.tetra_result = self.work_dir + '/tetra.summary.xls'
        self.tetra_pcoa = self.work_dir + '/tetra.pcoa.txt'

    def run_tetra(self):
        cmd = "{}checkm tetra -t 8 {} {}".format(self.checkm,self.option('input_genome').prop['path'],self.tetra_result)
        self.logger.info(cmd)
        self.logger.info("开始运行run_tetra")
        command = self.add_command("run_tetra", cmd)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.run_t_file(self.tetra_result,self.tetra_pcoa)
            self.logger.info("运行run_tetra完成")
        else:
            self.set_error("运行run_tetra运行出错!")

    def run_t_file(self,file,output):
        a = pd.read_table(file, sep="\t")
        b = a.T
        b.to_csv(output,sep="\t",header=0)

    def set_output(self):
        for i in ['tetra.summary.xls']:
            if os.path.exists(self.output_dir + "/" + i):
                os.remove(self.output_dir + "/" + i)
            os.link(self.work_dir + "/" + i,self.output_dir + "/" + i)
        self.option("out",self.tetra_pcoa)

    def run(self):
        """
        运行
        """
        super(CheckmTraTool, self).run()
        self.run_tetra()
        self.set_output()
        self.end()
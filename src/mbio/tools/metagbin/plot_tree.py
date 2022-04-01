#!/usr/bin/env python
# -*- coding: utf-8 -*-
# __author__ = 'gaohao'

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os,shutil,re
from mbio.packages.metagbin.common_function import link_dir,bin_rename

class PlotTreeAgent(Agent):
    """
    用于metagbin构建进化树，序列数量比较大
    version 1.0
    author: gaohao
    last_modify: 2019.03.04
    """

    def __init__(self, parent):
        super(PlotTreeAgent, self).__init__(parent)
        options = [
            {"name": "seq_fa", "type": "infile","format": "sequence.fasta"},  #
            {"name": "method", "type": "string", "default": "pro"},  # 构建进化树方法。蛋白 or 核苷酸
        ]
        self.add_option(options)

    def check_options(self):
        """
        检测参数是否正确
        """
        if not self.option('seq_fa').is_set:
            raise OptionError("组装的fasta文件不存在！")
        if self.option('method') not in ['pro','nuc']:
            raise OptionError("请提供准确的构建进化树的序列方法！")

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 4
        self._memory ='60G'

    def end(self):
        super(PlotTreeAgent, self).end()


class PlotTreeTool(Tool):
    """
    version 1.0
    """
    def __init__(self, config):
        super(PlotTreeTool, self).__init__(config)
        self.path =self.config.SOFTWARE_DIR + "/gcc/5.1.0/bin:" + self.config.SOFTWARE_DIR + "/bioinfo/align/mafft-7.299-with-extensions/bin/:" + self.config.SOFTWARE_DIR + "/bioinfo/phylogenetic/trimal-trimAl/source:"
        self.LD_LIBRARY_PATH = self.config.SOFTWARE_DIR + "/gcc/5.1.0/lib64"
        self.set_environ(PATH=self.path,LD_LIBRARY_PATH=self.LD_LIBRARY_PATH)
        self.fasta = self.option("seq_fa").prop['path']
        self.mafft = "../../../../../.." + self.config.PACKAGE_DIR + "/metagbin/mafft.sh"
        self.mafft_sof = self.config.SOFTWARE_DIR + "/bioinfo/align/mafft-7.299-with-extensions/bin/mafft-fftnsi"
        self.trimal = "/bioinfo/phylogenetic/trimal-trimAl/source/trimal"
        self.iqtree = "bioinfo/compare_genome/software/iqtree-1.6.12-Linux/bin/iqtree"
        self.aln = self.work_dir + "/mafft.aln.fas"
        self.trimal_fa = self.work_dir + "/trimal.aln.fas"

    def run_mafft(self):
        cmd = "{} {} {} {}".format(self.mafft,self.mafft_sof,self.fasta,self.aln)
        self.logger.info(cmd)
        self.logger.info("开始运行run_mafft")
        command = self.add_command("run_mafft", cmd)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("运行run_mafft完成")
        else:
            self.set_error("运行run_mafft运行出错!")

    def run_trimal(self):
        cmd = "{} -in {} -out {} -automated1".format(self.trimal,self.aln,self.trimal_fa)
        self.logger.info(cmd)
        self.logger.info("开始运行run_trimal")
        command = self.add_command("run_trimal", cmd)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("运行run_trimal完成")
        else:
            self.set_error("运行run_trimal运行出错!")

    def run_nuc_raxml(self):
        if os.path.exists(self.work_dir + "/temp"):
            shutil.rmtree(self.work_dir + "/temp")
        os.mkdir(self.work_dir + "/temp")
        cmd = "{} -s {} -pre {} -m HKY -nt AUTO -ntmax 8 -bb 1000 -mem 60G -st DNA".format(self.iqtree,self.trimal_fa,self.work_dir + "/temp/all")
        self.logger.info(cmd)
        self.logger.info("开始运行run_nuc_iqtree")
        command = self.add_command("run_nuc_iqtree", cmd)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("运行run_nuc_iqtree完成")
        else:
            self.set_error("运行run_nuc_iqtree运行出错!")

    def run_pro_raxml(self):
        if os.path.exists(self.work_dir + "/temp"):
            shutil.rmtree(self.work_dir + "/temp")
        os.mkdir(self.work_dir + "/temp")
        cmd = "{} -s {} -pre {} -m LG -nt AUTO -ntmax 8 -bb 1000 -mem 60G -st AA".format(self.iqtree,self.trimal_fa,self.work_dir + "/temp/all")
        self.logger.info(cmd)
        self.logger.info("开始运行run_pro_iqtree")
        command = self.add_command("run_pro_iqtree", cmd)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("运行run_pro_iqtree完成")
        else:
            self.set_error("运行run_pro_iqtree运行出错!")

    def set_output(self):
        if os.path.exists(self.output_dir + "/phylo_tree.nwk"):
            os.remove(self.output_dir + "/phylo_tree.nwk")
        os.link(self.work_dir + "/temp/all.treefile",self.output_dir + "/phylo_tree.nwk")

    def run(self):
        """
        运行
        """
        super(PlotTreeTool, self).run()
        if self.option("method") in ["nuc"]:
            self.run_mafft()
            self.run_trimal()
            self.run_nuc_raxml()
        elif self.option("method") in ["pro"]:
            self.run_mafft()
            self.run_trimal()
            self.run_pro_raxml()
        self.set_output()
        self.end()

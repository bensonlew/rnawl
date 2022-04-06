# -*- coding: utf-8 -*-


import os
from biocluster.core.exceptions import OptionError
from biocluster.agent import Agent
from biocluster.tool import Tool
import subprocess



class PhyloAlphaAgent(Agent):
    """
    """

    def __init__(self, parent):
        super(PhyloAlphaAgent, self).__init__(parent)
        options = [
            {"name": "table", "type": "infile","format":"tool_lab.simple"},
            {"name": "opt_file", "type": "infile","format":"tool_lab.simple"},
            {"name": "file_type", "type": "string","default":"is_tree"},
            {"name": "index", "type": "string","default":"pd"}
        ]
        self.add_option(options)

    def check_options(self):
        """
        检查参数
        :return:
        """
        if not self.option('table').is_set:
            raise OptionError('必须输入表格')
        if not self.option('opt_file').is_set:
            raise OptionError('必须输入进化树文件或序列文件')
        return True

    def set_resource(self):
        """
        设置所需资源，需在子类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 1
        self._memory = "10G"

    def end(self):
        super(PhyloAlphaAgent, self).end()


class PhyloAlphaTool(Tool):
    def __init__(self, config):
        super(PhyloAlphaTool, self).__init__(config)
        self.python_path =  'miniconda2/bin/python'
        self.r_path = self.config.SOFTWARE_DIR + '/program/miniconda3/envs/R/bin/Rscript'
        self.mafft =self.config.SOFTWARE_DIR + '/bioinfo/align/mafft-7.299-with-extensions/bin/mafft'
        self.fastree = self.config.SOFTWARE_DIR + '/bioinfo/phylogenetic/fasttree2.1.9/FastTreeMP'

    def run(self):
        """
        运行
        :return:
        """
        super(PhyloAlphaTool, self).run()
        if self.option('file_type') == 'is_tree':
            tree_path = self.option('opt_file').path
        else:
            fasta_file = self.option('opt_file').path
            cmd1 = self.mafft + ' --thread 10 %s > phylo.align'%fasta_file
            command = subprocess.Popen(cmd1, shell=True)
            command.communicate()
            if command.returncode == 0:
                self.logger.info("完成比对！")
            else:
                self.set_error("mafft运行出错！")

            cmd2 = self.fastree + ' -nt phylo.align > phylo.tre'
            command = subprocess.Popen(cmd2, shell=True)
            command.communicate()
            if command.returncode == 0:
                self.logger.info("完成进化！")
            else:
                self.set_error("进化运行出错！")

            tree_path = self.work_dir + '/phylo.tre'


        cmd =self.python_path + ' ' + self.config.PACKAGE_DIR + "/tool_lab/phylo_alpha.py %s %s %s %s "%(
            self.option('index'), tree_path, self.option('table').path, self.r_path)


        self.logger.info(cmd)
        self.logger.info("开始运行")

        command = self.add_command("run_phylo_alpha", cmd)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("运行run_phylo_alpha完成")
        else:
            self.set_error("运行run_phylo_alpha运行出错!")

        os.link(self.work_dir+'/%s.xls'%self.option('index'), self.output_dir+'/%s.xls'%self.option('index'))

        self.end()



# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
# last_modify : 2019.10.14

import re,os,shutil
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from mbio.packages.align.blast.xml2table import xml2table
from mbio.packages.bac_comp_genome.common_function import link_file


class CompareTreeAgent(Agent):
    """
    iqtree进行物种树和基因树的比较
    """

    def __init__(self, parent):
        super(CompareTreeAgent, self).__init__(parent)
        options = [
            {"name": "species_tree", "type": "infile", "format": "graph.newick_tree"},  # 物种树文件夹
            {"name": "gene_tree", "type": "infile", "format": "graph.newick_tree"}, #基因树文件夹
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("species_tree"):
            raise OptionError("必须输入species_tree文件！")
        if not self.option("gene_dir"):
            raise OptionError("必须输入gene_tree文件！")
        return True

    def set_resource(self):
        self._cpu = 2
        self._memory = '20G'

    def end(self):
        super(CompareTreeAgent, self).end()


class CompareTreeTool(Tool):
    def __init__(self, config):
        super(CompareTreeTool, self).__init__(config)
        self._version = "1.0"
        self.sof_path = self.config.SOFTWARE_DIR + "/bioinfo/compare_genome/software/iqtree-1.6.12-Linux/bin/iqtree"
        self.species = self.option("species_tree").prop['path']
        self.sh = "../../../../../.." + self.config.PACKAGE_DIR + "/bac_comp_genome/roary.sh"
        self.gene= self.option("gene_tree").prop['path']
        self.outfile = self.output_dir + "/compare_result.xls"

    def run(self):
        """
        运行
        :return:
        """
        super(CompareTreeTool, self).run()
        self.run_tree()
        self.set_output()
        self.end()

    def run_tree(self):
        cmd = '{} {} {} {} {}'.format(self.sh, self.sof_path, self.species, self.gene, self.outfile)
        command = self.add_command("run_tree", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_tree运行完成！")
        else:
            self.set_error("run_tree运行完成运行出错!")

    def set_output(self):
        pass

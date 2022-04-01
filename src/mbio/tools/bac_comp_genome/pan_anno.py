# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'
# version 1.0
# last_modify: 20119.10.10

import os
import re, shutil
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from mbio.packages.bac_comp_genome.common_function import link_dir

class PanAnnoAgent(Agent):
    """

    """
    def __init__(self, parent):
        super(PanAnnoAgent, self).__init__(parent)
        options = [
            {"name": "anno_file", "type": "infile", "format": "sequence.profile_table"},  # 物种注释总览表
            {"name": "cluster", "type": "infile", "format": "sequence.profile_table"},  # 聚类结果表
        ]
        self.add_option(options)
        self._memory_increase_step = 50  # 每次重运行增加内存50G by qingchen.zhang


    def check_options(self):
        if not self.option("cluster").is_set :
            raise OptionError("必须设置参数cluster文件!")
        if not self.option("anno_file").is_set:
            raise OptionError("必须设置参数anno_file文件!")

    def set_resource(self):
        self._cpu = 1
        number = os.path.getsize(self.option("anno_file").prop['path'])
        total_memory = number *2 / (1000*1000)
        if total_memory > 100:
            total_memory = 80
        self._memory = '{}G'.format(total_memory)

    def end(self):
        super(PanAnnoAgent, self).end()

class PanAnnoTool(Tool):
    def __init__(self, config):
        super(PanAnnoTool, self).__init__(config)
        self.python = "/program/Python/bin/python"
        self.merge_cluster = self.config.PACKAGE_DIR + "/bac_comp_genome/merge_cluster.py"
        self.out = self.work_dir + '/cluster_anno.xls'

    def run_anno(self):
        """
        将聚类表和注释表进行合并
        :return:
        """
        self.logger.info("开始对聚类和注释结果进行合并")
        if os.path.exists(self.out):
            os.remove(self.out)
        cmd = "{} {} -i {} -a {} -o {}".format(self.python, self.merge_cluster, self.option("cluster").prop['path'], self.option("anno_file").prop['path'], self.out)
        command = self.add_command("run_anno", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_anno运行完成！")
        else:
            self.set_error("run_anno运行出错!")

    def set_output(self):
        """
        设置结果文件目录
        :return:
        """
        self.logger.info("开始设置结果文件目录")
        if os.path.exists(os.path.join(self.output_dir, 'cluster_anno.xls')):
            os.remove(os.path.join(self.output_dir, 'cluster_anno.xls'))
        os.link(self.out, os.path.join(self.output_dir, 'cluster_anno.xls'))

    def run(self):
        super(PanAnnoTool, self).run()
        self.run_anno()
        self.set_output()
        self.end()

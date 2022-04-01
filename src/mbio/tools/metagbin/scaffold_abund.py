#!/usr/bin/env python
# -*- coding: utf-8 -*-
# __author__ = 'gaohao'

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os

class ScaffoldAbundAgent(Agent):
    """
    用于计算scaffolds的测序深度
    version 1.0
    author: gaohao
    last_modify: 2019.01.07
    """

    def __init__(self, parent):
        super(ScaffoldAbundAgent, self).__init__(parent)
        options = [
            {"name": "bam", "type": "string"},
            {"name": "num", "type": "int"},
        ]
        self.add_option(options)
        self.list =[]


    def check_options(self):
        """
        检测参数是否正确
        """
        if not self.option('bam'):
            raise OptionError("reads mapping的bam文件字段不存在！")


    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 2
        self._memory = '20G'

    def end(self):
        super(ScaffoldAbundAgent, self).end()


class ScaffoldAbundTool(Tool):
    """
    version 1.0
    """
    def __init__(self, config):
        super(ScaffoldAbundTool, self).__init__(config)
        self.bam =self.option('bam')
        self.num = self.option('num')
        self.summarize_depth ="/bioinfo/metaGenomic/metabat/jgi_summarize_bam_contig_depths"
        self.depth =self.work_dir + '/depth' + str(self.num) + '.txt'

    def run_summarize_depth(self):
        cmd = "{} --outputDepth {} {} ".format(self.summarize_depth, self.depth, self.bam)
        self.logger.info(cmd)
        self.logger.info("开始运行run_summarize_depth")
        command = self.add_command("run_summarize_depth", cmd)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("运行run_summarize_depth完成")
        else:
            self.set_error("运行run_summarize_depth运行出错!")

    def set_output(self):
        for i in ['depth' + str(self.num) + '.txt']:
            if os.path.exists(self.output_dir + "/" + i):
                os.remove(self.output_dir + "/" + i)
            os.link(self.work_dir + "/" + i,self.output_dir + "/" + i)

    def run(self):
        """
        运行
        """
        super(ScaffoldAbundTool, self).run()
        self.run_summarize_depth()
        self.set_output()
        self.end()
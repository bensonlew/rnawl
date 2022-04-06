# -*- coding: utf-8 -*-
# __author__ = 'haidong.gu'
# __last_modify__ = '2019/5/6'
from biocluster.agent import Agent
from biocluster.tool import Tool
import os
from biocluster.core.exceptions import OptionError


class ScafMeanCovAgent(Agent):
    """
    DemoAgent:
    version 1.0
    """

    def __init__(self, parent):
        super(ScafMeanCovAgent, self).__init__(parent)
        options = [
            {"name": "bam_file", "type": "infile", "format": "align.bwa.bam"},
            {"name": "depth", "type": "outfile", "format": "sequence.profile_table"}
        ]
        self.add_option(options)

    def check_options(self):
        """
        检查参数是否正确
        """
        # modified check_option
        return True

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 2
        self._memory = "10G"


class ScafMeanCovTool(Tool):
    def __init__(self, config):
        super(ScafMeanCovTool, self).__init__(config)
        self.script = self.config.PACKAGE_DIR + '/bacgenome/depth_mean.py'
        self.python_path = '/miniconda2/bin/python'
        self.depth = self.work_dir + "/" + 'cov.txt'
        self.depth_mean = self.work_dir + "/" + "cov_mean.txt"
        self.samtools = self.config.SOFTWARE_DIR + "/bioinfo/align/samtools-1.4/bin/samtools"

    def run_coverage(self):
        """
        description
        :return:
        """
        cmd = self.samtools + " depth " + self.option("bam_file").prop["path"] + " > " + self.depth
        command = self.add_command("coverage", cmd, shell=True).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("位点计算完成")
        else:
            self.set_error("位点计算出错")

    def run_mean(self):
        cmd = self.python_path + " " + self.script + " -n 100 -i " + self.depth + " -o " + self.depth_mean
        command = self.add_command("mean", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("mean计算完成")
        else:
            self.set_error("mean计算出错")

    def set_output(self):
        """
        设置输出文件路径
        :return:
        """
        os.link(self.depth_mean, self.output_dir + "/cov_mean.txt")
        self.option("depth").set_path(self.depth_mean)

    def run(self):
        super(ScafMeanCovTool, self).run()
        self.run_coverage()
        self.run_mean()
        self.set_output()
        self.end()
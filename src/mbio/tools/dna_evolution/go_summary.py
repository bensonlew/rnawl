# -*- coding: utf-8 -*-
# __author__ = 'wentianliu'
# last modify 20181102

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import datetime
import random
import os
import re


class GoSummaryAgent(Agent):
    """
    """
    def __init__(self, parent):
        super(GoSummaryAgent, self).__init__(parent)
        options = [
            {"name": "pop_summary", "type": "infile", "format": "dna_evolution.pop_summary"}
        ]
        self.add_option(options)
        self.step.add_steps('script')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.script.start()
        self.step.update()

    def step_end(self):
        self.step.script.finish()
        self.step.update()

    def check_options(self):
        if not self.option('pop_summary'):
            raise OptionError('必须输入:pop_summary')

    def set_resource(self):
        """
        运行所需资源
        """
        self._cpu = 2
        self._memory = '10G'

    def end(self):
        super(GoSummaryAgent, self).end()


class GoSummaryTool(Tool):
    def __init__(self, config):
        super(GoSummaryTool, self).__init__(config)
        self.perl_path = 'miniconda2/bin/perl'
        self.go_summary_path = self.config.PACKAGE_DIR + "/dna_evolution/go-summary.pl"
        self.obo_file = self.config.SOFTWARE_DIR + "/database/dna_wgs_geneome/obo_file/go-basic.obo"

    def run_go_summary(self):
        """
        obofile:库文件，每个项目都用这份。
        :return:
        """
        cmd = "{} {} -annotation {} -obofile {} -out {}"\
            .format(self.perl_path, self.go_summary_path, self.option('pop_summary').prop['path'],
                    self.obo_file, os.path.join(self.output_dir, "pop"))
        command = self.add_command("go_summary", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("go_summary运行成功")
        else:
            self.set_error("go_summary运行失败")

    def run(self):
        super(GoSummaryTool, self).run()
        self.run_go_summary()
        self.end()

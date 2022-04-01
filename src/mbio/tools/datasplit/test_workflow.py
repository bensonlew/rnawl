# -*- coding: utf-8 -*-
# __author__ = 'xuting'

"""bcl2fastq 工具 """
import os
import errno
import re
import xml.etree.ElementTree as ET
from biocluster.tool import Tool
from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError
from mbio.packages.datasplit.miseq_split import code2index


class TestWorkflowAgent(Agent):
    """
    bcl2fastq
    version 2.17
    """
    def __init__(self, parent=None):
        super(TestWorkflowAgent, self).__init__(parent)
        self._run_mode = "ssh1"
        options = [
            {'name': 'sample_info', 'type': "string"}
        ]
        self.add_option(options)
        self.step.add_steps("test")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.test.start()
        self.step.update()

    def stepfinish(self):
        self.step.test.finish()
        self.step.update()

    def check_options(self):
        """
        参数检测
        """
        return True

    def set_resource(self):
        """
        设置所需要的资源
        """
        self._cpu = 5
        self._memory = ''

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        super(TestWorkflowAgent, self).end()


class TestWorkflowTool(Tool):
    """
    """
    def __init__(self, config):
        super(TestWorkflowTool, self).__init__(config)

    def bcl2fastq(self):
        """
        运行bcl2fastq
        """
        os.link(self.option("sample_info") + '/RTAComplete.txt', self.output_dir + '/RTAComplete.txt')

    def run(self):
        super(TestWorkflowTool, self).run()
        self.bcl2fastq()
        self.end()

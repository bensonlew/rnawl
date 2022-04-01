# -*- coding: utf-8 -*-
# __author__ = 'xuxi'

import os
import unittest
from Bio import SeqIO
from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError
from biocluster.tool import Tool
import json


class DrawStaticplotAgent(Agent):
    """
    绘制peperr静态散点图
    """

    def __init__(self, parent):
        super(DrawStaticplotAgent, self).__init__(parent)
        options = [
            {"name": "inputfile", "type": "infile", "format": "ref_rna_v2.common"},
            {"name": "outdir", "type": "string", "default": None}

        ]
        self.add_option(options)

    def check_options(self):
        """
        重写参数检测函数
        :return:
        """
        if not self.option("inputfile").is_set:
            raise OptionError("必须提供inputfile输入文件")
        self.set_resource()
        return True

    def set_resource(self):
        """
        设置所需资源，需在类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 5
        self._memory = '10G'

    def end(self):
        super(DrawStaticplotAgent, self).end()


class DrawStaticplotTool(Tool):
    def __init__(self, config):
        super(DrawStaticplotTool, self).__init__(config)
        software_dir = self.config.SOFTWARE_DIR
        self.python3 = "/bioinfo/tool_lab/miniconda3/bin/python"
        self.plot_py = self.config.PACKAGE_DIR + '/itraq_and_tmt/peperr_image.py'


    def run(self):
        """
        运行
        :return:
        """
        super(DrawStaticplotTool, self).run()
        self.plot()
        self.end()

    def plot(self):
        cmd = '{} {} '.format(self.python3, self.plot_py)
        if self.option("outdir"):
            cmd += '{} {} '.format(self.option("inputfile").prop['path'], self.option("outdir"))
        else:
            cmd += '{} {} '.format(self.option("inputfile").prop['path'], self.output_dir)
        cmd_name = 'draw static plot'
        command = self.add_command(cmd_name, cmd, ignore_error=True)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("{} Finished successfully".format(cmd_name))
        else:
            self.logger.info("{} unfinished command".format(cmd_name))



class TestFunction(unittest.TestCase):
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "chart" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "itraq_and_tmt.chart",
            "instant": False,
            "options": dict(
                file_json="/mnt/ilustre/users/sanger-dev/sg-users/liubinxu/test_ref_rna_v2/chart_workflow.json"
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()

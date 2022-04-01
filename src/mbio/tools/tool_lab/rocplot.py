# -*- coding: utf-8 -*-
# __author__ = 'xuxi'

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import unittest
import re


class RocplotAgent(Agent):
    """
    """
    def __init__(self, parent):
        super(RocplotAgent, self).__init__(parent)
        options = [
            {"name": "infile", "type": "infile", "format": "ref_rna_v2.common"}
        ]
        self.add_option(options)
        self.step.add_steps('rocplot')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.rocplot.start()
        self.step.update()

    def step_end(self):
        self.step.rocplot.finish()
        self.step.update()

    def check_options(self):
        if not self.option('infile').is_set:
            raise OptionError('必须输入roc文件。')
        return True

    def set_resource(self):
        """
        运行所需资源
        """
        self._cpu = 1
        self._memory = "5G"

    def end(self):
        super(RocplotAgent, self).end()


class RocplotTool(Tool):
    def __init__(self, config):
        super(RocplotTool, self).__init__(config)
        self.python = 'bioinfo/tool_lab/miniconda3/bin/python'
        self.rocplot_py = os.path.join(self.config.PACKAGE_DIR, 'tool_lab/rocplot.py')

    def do_rocplot(self):
        """
        """
        cmd = "{} {} ".format(self.python, self.rocplot_py)
        cmd += '-infile {} '.format(self.option('infile').prop['path'])
        cmd += '-outdir {} '.format(self.output_dir)
        self.logger.info(cmd)
        self.logger.info("开始进行绘制rocplot")
        command = self.add_command("rocplot", cmd)
        command.run()  # 必须小写，
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("绘制rocplot完成")
        else:
            self.set_error("绘制rocplot出错！")

    def run(self):
        super(RocplotTool, self).run()
        self.do_rocplot()
        self.end()


class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run script to do test.
    """

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "rocplot_{}_{}".format(random.randint(1000, 9999), random.randint(1000, 9999)),
            "type": "tool",
            "name": "tool_lab.rocplot",
            "instant": False,
            "options": dict(
                infile='/mnt/lustre/users/sanger-dev/sg-users/xuxi/rocplot_test_dir/test_data.txt',
                target_lab='species',
                show_apart="heat-tree"
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == "__main__":
    suite = unittest.TestSuite()
    suite.addTests([TestFunction("test")])
    unittest.TextTestRunner(verbosity=2).run(suite)
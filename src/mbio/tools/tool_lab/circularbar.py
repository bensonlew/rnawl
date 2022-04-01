# -*- coding: utf-8 -*-
# __author__ = 'xuxi'

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import unittest
import re


class CircularbarAgent(Agent):
    """
    """
    def __init__(self, parent):
        super(CircularbarAgent, self).__init__(parent)
        options = [
            {"name": "infile", "type": "infile", "format": "ref_rna_v2.common"},# 第一列为id列，第二列为分组列（可有可无），最后一列为数值列
        ]
        self.add_option(options)
        self.step.add_steps('circularbar')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.circularbar.start()
        self.step.update()

    def step_end(self):
        self.step.circularbar.finish()
        self.step.update()

    def check_options(self):
        if not self.option('infile').is_set:
            raise OptionError('必须输入表达量文件。')
        return True

    def set_resource(self):
        """
        运行所需资源
        """
        self._cpu = 1
        self._memory = "5G"

    def end(self):
        super(CircularbarAgent, self).end()


class CircularbarTool(Tool):
    def __init__(self, config):
        super(CircularbarTool, self).__init__(config)
        # self.gcc = self.config.SOFTWARE_DIR + '/gcc/5.1.0/bin'
        # self.gcc_lib = self.config.SOFTWARE_DIR + '/gcc/5.1.0/lib64'
        # self.set_environ(PATH=self.gcc, LD_LIBRARY_PATH=self.gcc_lib)
        self.rscript = 'bioinfo/rna/miniconda2/bin/Rscript'
        self.circularbar_R = os.path.join(self.config.PACKAGE_DIR, 'tool_lab/circular_bar.r')

    def do_circularbar(self):
        """
        """
        cmd = "{} {} ".format(self.rscript, self.circularbar_R)
        cmd += '--infile {} '.format(self.option('infile').prop['path'])
        cmd += '--outdir {} '.format(self.output_dir)
        self.logger.info(cmd)
        self.logger.info("开始绘制circularbar")
        command = self.add_command("circularbar", cmd)
        command.run()  # 必须小写，
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("绘制circularbar完成")
        else:
            self.set_error("绘制circularbar出错！")

    def run(self):
        super(CircularbarTool, self).run()
        self.do_circularbar()
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
            "id": "circularbar_{}_{}".format(random.randint(1000, 9999), random.randint(1000, 9999)),
            "type": "tool",
            "name": "tool_lab.circularbar",
            "instant": False,
            "options": dict(
                infile='/mnt/lustre/users/sanger-dev/sg-users/xuxi/cicular_bar_test_dir/circular_bar_test_data.tsv',
                # infile='/mnt/lustre/users/sanger-dev/sg-users/xuxi/cicular_bar_test_dir/circular_bar_test_data_no_group.tsv',
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == "__main__":
    suite = unittest.TestSuite()
    suite.addTests([TestFunction("test")])
    unittest.TextTestRunner(verbosity=2).run(suite)
# -*- coding: utf-8 -*-
# __author__ = 'xuxi'

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import unittest
import re


class TreeheatAgent(Agent):
    """
    """
    def __init__(self, parent):
        super(TreeheatAgent, self).__init__(parent)
        options = [
            {"name": "infile", "type": "infile", "format": "ref_rna_v2.common"},
            {"name": "target_lab", "type": "string", 'default': ""},
            {"name": "show_apart", "type": "string", 'default': "heat-tree"}
        ]
        self.add_option(options)
        self.step.add_steps('treeheat')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.treeheat.start()
        self.step.update()

    def step_end(self):
        self.step.treeheat.finish()
        self.step.update()

    def check_options(self):
        if not self.option('infile').is_set:
            raise OptionError('必须输入表达量文件。')
        if not self.option('target_lab'):
            raise OptionError('必须输入目标列名。')
        return True

    def set_resource(self):
        """
        运行所需资源
        """
        self._cpu = 1
        self._memory = "5G"

    def end(self):
        super(TreeheatAgent, self).end()


class TreeheatTool(Tool):
    def __init__(self, config):
        super(TreeheatTool, self).__init__(config)
        # self.gcc = self.config.SOFTWARE_DIR + '/gcc/5.1.0/bin'
        # self.gcc_lib = self.config.SOFTWARE_DIR + '/gcc/5.1.0/lib64'
        # self.set_environ(PATH=self.gcc, LD_LIBRARY_PATH=self.gcc_lib)
        self.rscript = 'bioinfo/ref_rna_v3/clusterprofile_4.1/miniconda3/bin/Rscript'
        self.treeheat_R = os.path.join(self.config.PACKAGE_DIR, 'tool_lab/treeheat.r')

    def do_treeheat(self):
        """
        """
        cmd = "{} {} ".format(self.rscript, self.treeheat_R)
        cmd += '--infile {} '.format(self.option('infile').prop['path'])
        cmd += '--outdir {} '.format(self.output_dir)
        cmd += '--target_lab {} '.format(self.option('target_lab'))
        cmd += '--show_apart {} '.format(self.option('show_apart'))
        self.logger.info(cmd)
        self.logger.info("开始进行绘制热图决策树treeheat")
        command = self.add_command("treeheat", cmd)
        command.run()  # 必须小写，
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("热图决策树treeheat完成")
        else:
            self.set_error("热图决策树treeheat出错！")

    def run(self):
        super(TreeheatTool, self).run()
        self.do_treeheat()
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
            "id": "treeheat_{}_{}".format(random.randint(1000, 9999), random.randint(1000, 9999)),
            "type": "tool",
            "name": "tool_lab.treeheat",
            "instant": False,
            "options": dict(
                infile='/mnt/lustre/users/sanger-dev/sg-users/xuxi/treeheat_test_dir/test_data.txt',
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
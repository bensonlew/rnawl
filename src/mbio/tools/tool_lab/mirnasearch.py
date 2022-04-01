# -*- coding: utf-8 -*-
# __author__ = 'xuxi'

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import unittest
import re


class MirnasearchAgent(Agent):
    """
    """
    def __init__(self, parent):
        super(MirnasearchAgent, self).__init__(parent)
        options = [
            {"name": "speice", "type": "string", "format": "Human"},
            {"name": "gene", "type": "string", 'default': ""},
            {"name": "mirna", "type": "string", 'default': ""},
            {"name": "database", "type": "string", 'default': "all"}
        ]
        self.add_option(options)
        self.step.add_steps('mirnasearch')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.mirnasearch.start()
        self.step.update()

    def step_end(self):
        self.step.mirnasearch.finish()
        self.step.update()

    def check_options(self):
        if not self.option('speice'):
            raise OptionError('必须搜索物种名。')
        return True

    def set_resource(self):
        """
        运行所需资源
        """
        self._cpu = 5
        self._memory = "5G"

    def end(self):
        super(MirnasearchAgent, self).end()


class MirnasearchTool(Tool):
    def __init__(self, config):
        super(MirnasearchTool, self).__init__(config)
        # self.gcc = self.config.SOFTWARE_DIR + '/gcc/5.1.0/bin'
        # self.gcc_lib = self.config.SOFTWARE_DIR + '/gcc/5.1.0/lib64'
        # self.set_environ(PATH=self.gcc, LD_LIBRARY_PATH=self.gcc_lib)
        self.python = 'program/miniconda3/bin/python3'
        self.miRTarBase_file = self.config.SOFTWARE_DIR + "/database/Tool_lab/mirnasearch/miRTarBase_MTI.txt"
        self.idmappdb = self.config.SOFTWARE_DIR + "/database/Tool_lab/mirnasearch/human_mouse_rat_id_db.txt"
        self.mirnasearch_py = os.path.join(self.config.PACKAGE_DIR, 'tool_lab/mirnasearch.py')

    def do_mirnasearch(self):
        """
        """
        cmd = "{} {} ".format(self.python, self.mirnasearch_py)
        cmd += '-mirtarbase {} '.format(self.miRTarBase_file)
        cmd += '-idmappdb {} '.format(self.idmappdb)
        cmd += '-speice {} '.format(self.option('speice'))
        if self.option('mirna'):
            cmd += '-mirna {} '.format(self.option('mirna'))
        if self.option('gene'):
            cmd += '-gene {} '.format(self.option('gene'))
        cmd += '-outfile {} '.format(os.path.join(self.output_dir,"result.txt"))
        cmd += '-database {} '.format(self.option('database'))
        self.logger.info(cmd)
        self.logger.info("开始进行mirnasearch")
        command = self.add_command("mirnasearch", cmd)
        command.run()  # 必须小写，
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("mirnasearch完成")
        else:
            self.set_error("mirnasearch出错！")

    def run(self):
        super(MirnasearchTool, self).run()
        self.do_mirnasearch()
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
            "id": "mirnasearch_{}_{}".format(random.randint(1000, 9999), random.randint(1000, 9999)),
            "type": "tool",
            "name": "tool_lab.mirnasearch",
            "instant": False,
            "options": dict(
                infile='/mnt/lustre/users/sanger-dev/sg-users/xuxi/mirnasearch_test_dir/test_data.txt',
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
# -*- coding: utf-8 -*-
# __author__ = 'zhangyitong'

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import unittest
import re


class ManhattanAgent(Agent):
    """
    """
    def __init__(self, parent):
        super(ManhattanAgent, self).__init__(parent)
        options = [
            {"name": "model_file", "type": "infile", "format": "ref_rna_v2.common"},
            {"name": "method", "type": "string", 'default': 'variants_index'},
        ]
        self.add_option(options)
        self.step.add_steps('manhattan_draw')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.manhattan_draw.start()
        self.step.update()

    def step_end(self):
        self.step.manhattan_draw.finish()
        self.step.update()

    def check_options(self):
        if not self.option('model_file').is_set:
            raise OptionError('必须设置输入模型结果文件。')
        return True

    def set_resource(self):
        """
        运行所需资源
        """
        self._cpu = 1
        self._memory = "20G"

    def end(self):
        super(ManhattanAgent, self).end()


class ManhattanTool(Tool):
    def __init__(self, config):
        super(ManhattanTool, self).__init__(config)
        self.gcc = self.config.SOFTWARE_DIR + '/gcc/5.1.0/bin'
        self.gcc_lib = self.config.SOFTWARE_DIR + '/gcc/5.1.0/lib64'
        self.set_environ(PATH=self.gcc, LD_LIBRARY_PATH=self.gcc_lib)
        self.program = {
            'rscript': 'bioinfo/rna/miniconda2/bin/Rscript',
            'perl': 'program/perl/perls/perl-5.24.0/bin/perl',
        }
        self.script = {
            'index': os.path.join(self.config.PACKAGE_DIR, 'tool_lab/bsa/manhattan-index.r'),
            'gvalue': os.path.join(self.config.PACKAGE_DIR, 'tool_lab/bsa/manhattan-G.r'),
        }
        self.file = {
            'out': os.path.join(self.output_dir, 'pop'),
        }

    def draw_index(self):
        """
        """
        cmd = "{} {} ".format(self.program['rscript'], self.script['index'])
        cmd += '--result {} '.format(self.option('model_file').prop['path'])
        cmd += '--output {}'.format(self.file['out'])
        self.logger.info(cmd)
        self.logger.info("开始进行绘制Manhattan图")
        command = self.add_command("draw_manhattan", cmd)
        command.run()  # 必须小写，
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("绘制Manhattan图完成！")
        else:
            self.set_error("绘制Manhattan图出错！")

    def draw_gvalue(self):
        """
        """
        cmd = "{} {} ".format(self.program['rscript'], self.script['gvalue'])
        cmd += '--result {} '.format(self.option('model_file').prop['path'])
        cmd += '--output {}'.format(self.file['out'])
        self.logger.info(cmd)
        self.logger.info("开始进行绘制Manhattan图")
        command = self.add_command("draw_manhattan", cmd)
        command.run()  # 必须小写，
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("绘制Manhattan图完成！")
        else:
            self.set_error("绘制Manhattan图出错！")

    def run(self):
        super(ManhattanTool, self).run()
        if self.option('method') == 'variants_index':
            self.draw_index()
        else:
            self.draw_gvalue()
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
            "id": "bsa_manhattan_{}_{}".format(random.randint(1000, 9999), random.randint(1000, 9999)),
            "type": "tool",
            "name": "tool_lab.bsa.manhattan",
            "instant": False,
            "options": dict(
                model_file='/mnt/ilustre/users/sanger-dev/sg-users/zhangyitong/BSA/demo/Index/pop.bootstrap.result',
                method='variants_index'
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == "__main__":
    suite = unittest.TestSuite()
    suite.addTests([TestFunction("test")])
    unittest.TextTestRunner(verbosity=2).run(suite)
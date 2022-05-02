# -*- coding: utf-8 -*-
# __author__ = 'xuxi'

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import unittest
import re


class CorrheatmapAgent(Agent):
    """
    """
    def __init__(self, parent):
        super(CorrheatmapAgent, self).__init__(parent)
        options = [
            {"name": "infile", "type": "infile", "format": "ref_rna_v2.common"},
            {"name": "type", "type": "string", 'default': 'full'},
            {"name": "method", "type": "string", 'default': "circle"},
            {"name": "color", "type": "string", 'default': "Spectral"},
            {"name": "corr_method", "type": "string", 'default': "pearson"},
            {"name": "corr_adjust", "type": "string", 'default': 'holm'},
            {"name": "show_significance", "type": "string", 'default': 'yes'},
            {"name": "in_significance", "type": "string", 'default': 'pch'},
            {"name": "show_coefficient", "type": "string", 'default': 'yes'}
        ]
        self.add_option(options)
        self.step.add_steps('corrheatmap')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.corrheatmap.start()
        self.step.update()

    def step_end(self):
        self.step.corrheatmap.finish()
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
        super(CorrheatmapAgent, self).end()


class CorrheatmapTool(Tool):
    def __init__(self, config):
        super(CorrheatmapTool, self).__init__(config)
        # self.gcc = self.config.SOFTWARE_DIR + '/gcc/5.1.0/bin'
        # self.gcc_lib = self.config.SOFTWARE_DIR + '/gcc/5.1.0/lib64'
        # self.set_environ(PATH=self.gcc, LD_LIBRARY_PATH=self.gcc_lib)
        self.rscript = 'miniconda2/bin/Rscript'
        self.corrheatmap_R = os.path.join(self.config.PACKAGE_DIR, 'tool_lab/corrheatmap.r')

    def do_corrheatmap(self):
        """
        """
        cmd = "{} {} ".format(self.rscript, self.corrheatmap_R)
        cmd += '--infile {} '.format(self.option('infile').prop['path'])
        cmd += '--outdir {} '.format(self.output_dir)
        cmd += '--type {} '.format(self.option('type'))
        cmd += '--method {} '.format(self.option('method'))
        cmd += '--color {} '.format(self.option('color'))
        cmd += '--corr_method {} '.format(self.option('corr_method'))
        cmd += '--corr_adjust {} '.format(self.option('corr_adjust'))
        cmd += '--show_significance {} '.format(self.option('show_significance'))
        cmd += '--in_significance {} '.format(self.option('in_significance'))
        cmd += '--show_coefficient {} '.format(self.option('show_coefficient'))
        self.logger.info(cmd)
        self.logger.info("开始进行绘制相关性热图corrheatmap")
        command = self.add_command("corrheatmap", cmd)
        command.run()  # 必须小写，
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("相关性热图corrheatmap完成")
        else:
            self.set_error("相关性热图corrheatmap出错！")

    def run(self):
        super(CorrheatmapTool, self).run()
        self.do_corrheatmap()
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
            "id": "corrheatmap_{}_{}".format(random.randint(1000, 9999), random.randint(1000, 9999)),
            "type": "tool",
            "name": "tool_lab.corrheatmap",
            "instant": False,
            "options": dict(
                infile='/mnt/ilustre/users/sanger-dev/sg-users/xuxi/corrheatmap_test_dir/test.tsv',
                type='full',
                method="circle",
                color="Spectral",
                corr_method="pearson",
                corr_adjust="holm",
                show_significance="yes",
                in_significance="pch",
                show_coefficient="yes"
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == "__main__":
    suite = unittest.TestSuite()
    suite.addTests([TestFunction("test")])
    unittest.TextTestRunner(verbosity=2).run(suite)
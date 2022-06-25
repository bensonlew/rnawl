# -*- coding: utf-8 -*-
# __author__ = 'zhangyitong'

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import unittest
import re


class RegionCalcAgent(Agent):
    """
    """
    def __init__(self, parent):
        super(RegionCalcAgent, self).__init__(parent)
        options = [
            {"name": "model_file", "type": "infile", "format": "ref_rna_v2.common"},
            {"name": "method", "type": "string", 'default': 'variants_index'},
        ]
        self.add_option(options)
        self.step.add_steps('region_calc')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.region_calc.start()
        self.step.update()

    def step_end(self):
        self.step.region_calc.finish()
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
        super(RegionCalcAgent, self).end()


class RegionCalcTool(Tool):
    def __init__(self, config):
        super(RegionCalcTool, self).__init__(config)
        self.program = {
            'rscript': 'program/R-3.3.1/bin/Rscript',
            'perl': 'miniconda2/bin/perl',
        }
        self.script = {
            'index': os.path.join(self.config.PACKAGE_DIR, 'tool_lab/bsa/region-index.pl'),
            'gvalue': os.path.join(self.config.PACKAGE_DIR, 'tool_lab/bsa/region-G.pl'),
        }
        self.file = {
            'out': os.path.join(self.output_dir, 'pop.result'),
        }

    def index_region(self):
        """
        """
        bootstrap = os.path.splitext(self.option('model_file').prop['path'])[0] + '.detail'
        cmd = "{} {} ".format(self.program['perl'], self.script['index'])
        cmd += '-input {} '.format(self.option('model_file').prop['path'])
        cmd += '-bootstrap {} '.format(bootstrap)
        cmd += '-output {}'.format(self.file['out'])
        self.logger.info(cmd)
        self.logger.info("开始进行关联位点统计")
        command = self.add_command("calc_region", cmd)
        command.run()  # 必须小写，
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("关联位点统计完成！")
        else:
            self.set_error("关联位点统计出错！")

    def gvalue_region(self):
        """
        """
        cmd = "{} {} ".format(self.program['perl'], self.script['gvalue'])
        cmd += '-pvalue 0.05 '
        cmd += '-input {} '.format(self.option('model_file').prop['path'])
        cmd += '-output {}'.format(self.file['out'])
        self.logger.info(cmd)
        self.logger.info("开始进行关联位点统计")
        command = self.add_command("calc_region", cmd)
        command.run()  # 必须小写，
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("关联位点统计完成！")
        else:
            self.set_error("关联位点统计出错！")

    def run(self):
        super(RegionCalcTool, self).run()
        if self.option('method') == 'variants_index':
            self.index_region()
        else:
            self.gvalue_region()
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
            "id": "bsa_region_calc_{}_{}".format(random.randint(1000, 9999), random.randint(1000, 9999)),
            "type": "tool",
            "name": "tool_lab.bsa.region_calc",
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
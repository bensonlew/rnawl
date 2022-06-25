# -*- coding: utf-8 -*-
# __author__ = 'zhangyitong'

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import unittest
import re


class ModelConstructionAgent(Agent):
    """
    """
    def __init__(self, parent):
        super(ModelConstructionAgent, self).__init__(parent)
        options = [
            {"name": "pop_index", "type": "infile", "format": "ref_rna_v2.common"},
            {"name": "method", "type": "string", 'default': 'variants_index'},
            {'name': 'win_size', 'type': 'int', 'default': 2000000},
            {'name': 'win_step', 'type': 'int', 'default': 100000},
            {"name": "model_file", "type": "outfile", "format": "ref_rna_v2.common"},
        ]
        self.add_option(options)
        self.step.add_steps('bsa_calc')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.bsa_calc.start()
        self.step.update()

    def step_end(self):
        self.step.bsa_calc.finish()
        self.step.update()

    def check_options(self):
        if not self.option("pop_index").is_set:
            raise OptionError("请设置pop_index")
        if not self.option('win_size') > 0:
            raise OptionError('必须设置输入大于0的window size。')
        if not self.option('win_step') > 0:
            raise OptionError('必须设置输入大于0的window step。')
        if not self.option('win_size') > self.option('win_step'):
            raise OptionError('window size必须大于window step。')
        return True

    def set_resource(self):
        """
        运行所需资源
        """
        self._cpu = 1
        self._memory = "20G"

    def end(self):
        super(ModelConstructionAgent, self).end()


class ModelConstructionTool(Tool):
    def __init__(self, config):
        super(ModelConstructionTool, self).__init__(config)
        self.program = {
            'rscript': 'program/R-3.3.1/bin/Rscript',
            'perl': 'miniconda2/bin/perl',
        }
        self.script = {
            'win_index': os.path.join(self.config.PACKAGE_DIR, 'tool_lab/bsa/slidingwin-index.r'),
            'bootstrap': os.path.join(self.config.PACKAGE_DIR, 'tool_lab/bsa/bootstrap-index.r'),
            'g_calc': os.path.join(self.config.PACKAGE_DIR, 'tool_lab/bsa/G-calc.r')
        }
        self.file = {
            'index_out': os.path.join(self.output_dir, 'pop'),
            'index_result': os.path.join(self.output_dir, 'pop.bootstrap.result'),
            'index_detail': os.path.join(self.output_dir, 'pop.bootstrap.detail'),
            'g_result': os.path.join(self.output_dir, 'pop.Gprime')
        }

    def run_win_index(self):
        """
        """
        cmd = "{} {} ".format(self.program['rscript'], self.script['win_index'])
        cmd += '--infile {} '.format(self.option('pop_index').prop['path'])
        cmd += '--outfile {} '.format(self.file['index_out'])
        cmd += '--winsize {} --stepsize {} '.format(self.option('win_size'), self.option('win_step'))
        cmd += '--method bp'
        self.logger.info(cmd)
        self.logger.info("开始进行sliding window index计算")
        command = self.add_command("sliding_window_index", cmd)
        command.run()  # 必须小写，
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("sliding window index计算完成！")
        else:
            self.set_error("sliding window index计算出错！")
        self.run_bootstrap()

    def run_bootstrap(self):
        """
        """
        cmd = "{} {} ".format(self.program['rscript'], self.script['bootstrap'])
        cmd += '--infile {} '.format(self.file['index_out'] + '.sliding.detail')
        cmd += '--bulksize 30 --outfile {} '.format(self.file['index_out'])
        cmd += '--bootstrap 1000 --popstruc F2 --qvalue 0.999'
        self.logger.info(cmd)
        self.logger.info("开始进行bootstrap计算")
        command = self.add_command("run_bootstrap", cmd).run()  # 必须小写，
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("bootstrap计算完成！")
        else:
            self.set_error("bootstrap计算出错！")
        if os.path.exists(self.file['index_result']):
            self.logger.info("设置model_file成功")
            self.option("model_file", self.file['index_result'])

    def run_gvalue(self):
        """
        """
        cmd = "{} {} ".format(self.program['rscript'], self.script['g_calc'])
        cmd += '--infile {} '.format(self.option('pop_index').prop['path'])
        cmd += '--outfile {} '.format(self.file['g_result'])
        cmd += '--winsize {} '.format(self.option('win_size'))
        cmd += '--method bp'
        self.logger.info(cmd)
        self.logger.info("开始进行G value计算")
        command = self.add_command("run_gvalue", cmd)
        command.run()  # 必须小写，
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("G value计算完成！")
        else:
            self.set_error("G value计算出错！")
        if os.path.exists(self.file['g_result']):
            self.logger.info("设置model_file成功")
            self.option("model_file", self.file['g_result'])

    def run(self):
        super(ModelConstructionTool, self).run()
        if self.option('method') == 'variants_index':
            self.run_win_index()
        else:
            self.run_gvalue()
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
            "id": "selective_sweep_filter_{}_{}".format(random.randint(1000, 9999), random.randint(1000, 9999)),
            "type": "tool",
            "name": "tool_lab.selective_sweep.vcftools_filter",
            "instant": False,
            "options": dict(
                vcf_path='/mnt/ilustre/users/sanger-dev/sg-users/zhangyitong/test/tool_052021/selective_sweep/final.vcf',
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == "__main__":
    suite = unittest.TestSuite()
    suite.addTests([TestFunction("test")])
    unittest.TextTestRunner(verbosity=2).run(suite)
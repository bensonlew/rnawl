# -*- coding: utf-8 -*-
# __author__ = 'zhangyitong'

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import pandas as pd
import os
import unittest
import glob
import re


class GwasAgent(Agent):
    """
    """
    def __init__(self, parent):
        super(GwasAgent, self).__init__(parent)
        options = [
            {"name": "plink", "type": "infile", "format": "ref_rna_v2.common_dir"},
            {"name": 'trait', 'type': 'infile', 'format': "ref_rna_v2.common"},
            {'name': 'chrom_map', "type": 'infile', 'format': 'ref_rna_v2.common'},
            {"name": "method", "type": "string", 'default': 'mlm'},
            {"name": 'alpha', "type": "float", 'default': 0.05},
        ]
        self.add_option(options)
        self.step.add_steps('gwas')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.gwas.start()
        self.step.update()

    def step_end(self):
        self.step.gwas.finish()
        self.step.update()

    def check_options(self):
        if not self.option('plink').is_set:
            raise OptionError('必须设置输入plink结果文件。')
        if not self.option('plink').is_set:
            raise OptionError('必须设置输入选择性状文件。')
        return True

    def set_resource(self):
        """
        运行所需资源
        """
        self._cpu = 2
        self._memory = "20G"

    def end(self):
        super(GwasAgent, self).end()


class GwasTool(Tool):
    def __init__(self, config):
        super(GwasTool, self).__init__(config)
        self.gcc = self.config.SOFTWARE_DIR + '/gcc/5.1.0/bin'
        self.gcc_lib = self.config.SOFTWARE_DIR + '/gcc/5.1.0/lib64'
        self.set_environ(PATH=self.gcc, LD_LIBRARY_PATH=self.gcc_lib)
        self.program = {
            'rscript': 'bioinfo/rna/miniconda2/bin/Rscript',
        }
        self.script = {
            'gwas': os.path.join(self.config.PACKAGE_DIR, 'tool_lab/gwas.r'),
        }
        self.file = {
            'out': os.path.join(self.output_dir, 'gwas_result.xls'),
        }

    def run_gwas(self):
        """
        """
        cmd = "{} {} ".format(self.program['rscript'], self.script['gwas'])
        cmd += '--plink {} '.format(os.path.join(self.option('plink').prop['path'], 'pop'))
        cmd += '--trait {} '.format(self.option('trait').prop['path'])
        cmd += '--method {} '.format(self.option('method'))
        cmd += '--output {} '.format(self.output_dir)
        cmd += '--threshold {} '.format(self.option('alpha'))
        cmd += '--chrom {}'. format(self.option('chrom_map').prop['path'])
        self.logger.info(cmd)
        self.logger.info("开始进行gwas计算")
        command = self.add_command("run_gwas", cmd)
        command.run()  # 必须小写，
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("gwas计算完成！")
        else:
            self.set_error("gwas计算出错！")

    def run(self):
        super(GwasTool, self).run()
        self.run_gwas()
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
            "id": "gwas_gwas_{}_{}".format(random.randint(1000, 9999), random.randint(1000, 9999)),
            "type": "tool",
            "name": "tool_lab.gwas.gwas",
            "instant": False,
            "options": dict(
                plink='/mnt/ilustre/users/sanger-dev/workspace/20210519/Single_gwas_plink_6872_6589/Plink/output',
                trait='/mnt/ilustre/users/sanger-dev/sg-users/zhangyitong/test/tool_052021/gwas/trait.txt'
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == "__main__":
    suite = unittest.TestSuite()
    suite.addTests([TestFunction("test")])
    unittest.TextTestRunner(verbosity=2).run(suite)
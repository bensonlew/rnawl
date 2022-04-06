# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import os
import unittest

import pandas as pd
from biocluster.agent import Agent
from biocluster.tool import Tool

from mbio.packages.whole_transcriptome.utils import runcmd


class UniformAgent(Agent):
    """
    last_modify: 2019.10.12
    """

    def __init__(self, parent):
        super(UniformAgent, self).__init__(parent)
        PROGRAM = ('DESeq2', 'edgeR', 'DEGseq', 'Limma', 'NOISeq', 'svaseqlimma')
        STAT_TYPE = ('pvalue', 'padjust')
        options = [
            {'name': 'method', 'type': 'string', 'default': PROGRAM[0]},
            {'name': 'input_dir', 'type': 'infile', 'format': 'whole_transcriptome.common_dir'},
            {'name': 'pvalue_padjust', 'type': 'string', 'default': STAT_TYPE[0]}
        ]
        self.add_option(options)

    def check_options(self):
        for k, v in self._options.items():
            self.logger.debug('{} -> {}'.format(k, v.value))
        else:
            return True

    def set_resource(self):
        self._cpu = 1
        self._memory = '10G'

    def end(self):
        super(UniformAgent, self).end()


class UniformTool(Tool):
    def __init__(self, config):
        super(UniformTool, self).__init__(config)
        self.program = {
            'python': 'miniconda2/bin/python',
        }
        self.script = {
            'uniform': os.path.join(self.config.PACKAGE_DIR, 'denovo_rna_v3/batch/uniform.py')
        }

    def run(self):
        super(UniformTool, self).run()
        self.run_uniform()
        self.set_output()
        self.end()

    def run_uniform(self):
        cmd = '{} {}'.format(self.program['python'], self.script['uniform'])
        cmd += ' -m {}'.format(self.option('method'))
        cmd += ' -i {}'.format(self.option('input_dir').path)
        cmd += ' -o {}'.format(self.output_dir)
        cmd += ' -p {}'.format(self.option('pvalue_padjust'))
        cmd_name = 'uniform'
        command = self.add_command(cmd_name, cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("{} 运行成功".format(cmd_name))
        elif command.return_code is None:
            self.logger.warn("运行{}出错，返回值为None，尝试重新运行".format(cmd_name))
            command.rerun()
            self.wait()
            if command.return_code is 0:
                self.logger.info("{} 运行成功".format(cmd_name))
            else:
                self.set_error("运行%s>>>%s出错", variables=(cmd_name, cmd), code="33704306")
        else:
            self.set_error("运行%s>>>%s出错", variables=(cmd_name, cmd), code="33704307")

    def set_output(self):
        pass


class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'uniform_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'tool',
            'name': 'denovo_rna_v3.batch.uniform',
            'instant': False,
            'options': {
                'input_dir': '/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/uniform/output',
                'method': 'DESeq2',
                'pvalue_padjust': 'padjust'
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)

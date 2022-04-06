# -*- coding: utf-8 -*-
# __author__ = 'zoujiaxun'

import os
import unittest
import pandas as pd

from biocluster.agent import Agent
from biocluster.tool import Tool
from mbio.packages.whole_transcriptome.utils import runcmd
from biocluster.core.exceptions import OptionError


class NomogramAgent(Agent):
    '''
    last_modify: 2020.07.23
    '''

    def __init__(self, parent):
        super(NomogramAgent, self).__init__(parent)
        options = [
            {'name': 'surv_table', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'factor_list', 'type': 'string'},
            {'name': 'time_list', 'type': 'string'},
            {'name': 'out_table', 'type': 'outfile', 'format': 'whole_transcriptome.common'},
        ]
        self.add_option(options)

    def check_options(self):
        pass

    def set_resource(self):
        self._cpu = 1
        self._memory = '10G'

    def end(self):
        super(NomogramAgent, self).end()


class NomogramTool(Tool):
    def __init__(self, config):
        super(NomogramTool, self).__init__(config)
        self.program = {
            'python': 'miniconda2/bin/python',
        }
        self.script = {
            'nomogram': os.path.join(self.config.PACKAGE_DIR, 'tool_lab/nomogram.py')
        }
        self.file = {
            'out_table': os.path.join(self.output_dir, "nomogram.txt")
        }

    def run(self):
        super(NomogramTool, self).run()
        self.run_nomogram()
        self.set_output()
        self.end()

    def run_nomogram(self):
        cmd = '{} {} -i {} -fl {} -tl {} -o {}'.format(self.program['python'], self.script['nomogram'],
                                                       self.option('surv_table').path, self.option('factor_list'),
                                                       self.option('time_list'), self.file['out_table'])
        cmd_name = 'run_nomogram'
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
                self.set_error("运行%s>>>%s出错", variables = (cmd_name, cmd), code = "33704306")
        else:
            self.set_error("运行%s>>>%s出错", variables = (cmd_name, cmd), code = "33704307")

    def set_output(self):
        # self.option('out_table').set_path(self.file['out_table'])
        pass


class TestFunction(unittest.TestCase):
    '''
    This is test for the tool. Just run this script to do test.
    '''

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'nomogram_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'tool',
            'name': 'tool_lab.nomogram',
            'instant': False,
            'options': {
                'surv_table': '/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/tool_lab/nomogram/lung.txt',
                'factor_list': 'age;sex;ph.karno',
                'time_list': '365;730'
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)



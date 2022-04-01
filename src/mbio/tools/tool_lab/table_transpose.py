# -*- coding: utf-8 -*-
# __author__ = 'zoujiaxun,qinjincheng'

import os
import unittest
import pandas as pd

from biocluster.agent import Agent
from biocluster.tool import Tool
from mbio.packages.whole_transcriptome.utils import runcmd
from biocluster.core.exceptions import OptionError


class TableTransposeAgent(Agent):
    '''
    last_modify: 2019.12.12
    '''

    def __init__(self, parent):
        super(TableTransposeAgent, self).__init__(parent)
        options = [
            {'name': 'table', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'sep', 'type': 'string'},
            {'name': 'transposed_table', 'type': 'outfile', 'format': 'whole_transcriptome.common'}

        ]
        self.add_option(options)

    def check_options(self):
        pass

    def set_resource(self):
        self._cpu = 10
        self._memory = '10G'

    def end(self):
        super(TableTransposeAgent, self).end()


class TableTransposeTool(Tool):
    def __init__(self, config):
        super(TableTransposeTool, self).__init__(config)
        self.program = {
            'python': 'program/Python/bin/python'
        }
        self.script = {
            'table_kit': os.path.join(self.config.PACKAGE_DIR, 'tool_lab/table_kit.py')
        }
        self.file = {
            'transposed_table': os.path.join(self.output_dir, '{}'.format(os.path.basename(self.option('table').path))),
        }

    def run(self):
        super(TableTransposeTool, self).run()

        self.run_transposed()
        self.set_output()
        self.end()

    def run_transposed(self):
        cmd = '{} {}'.format(self.program['python'], self.script['table_kit'])
        cmd += ' {}'.format('transpose')
        cmd += ' --input {}'.format(self.option('table').path)
        cmd += ' --sep {}'.format(self.option('sep'))
        cmd += ' --output {}'.format(self.file['transposed_table'])
        cmd_name = 'run_transposed'
        runcmd(self, cmd_name, cmd)


    def set_output(self):
        self.option('transposed_table').set_path(self.file['transposed_table'])

class TestFunction(unittest.TestCase):
    '''
    This is test for the tool. Just run this script to do test.
    '''

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'table_transpose_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'tool',
            'name': 'tool_lab.table_transpose',
            'instant': False,
            'options': {
                'table': '/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/venn.txt',
                'sep': 'tab',


            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)



# -*- coding: utf-8 -*-
# __author__ = 'zoujiaxun,qinjincheng'

import os
import unittest
import pandas as pd

from biocluster.agent import Agent
from biocluster.tool import Tool
from mbio.packages.whole_transcriptome.utils import runcmd
from biocluster.core.exceptions import OptionError


class TableStandardAgent(Agent):
    '''
    last_modify: 2019.12.12
    '''

    def __init__(self, parent):
        super(TableStandardAgent, self).__init__(parent)
        options = [
            {'name': 'table', 'type': 'infile', 'format': 'ref_rna_v2.common'},
            # {'name': 'manipulating', 'type': 'string'},
            {'name': 'sep', 'type': 'string'},
            {'name': 'observe', 'type': 'string'},
            {'name': 'guard', 'type': 'string'},
            {'name': 'ref_feature', 'type': 'string'},
            {'name': 'feature', 'type': 'string'},
            {'name': 'transposed_table', 'type': 'outfile', 'format': 'ref_rna_v2.common'}

        ]
        self.add_option(options)

    def check_options(self):
        pass

    def set_resource(self):
        self._cpu = 10
        self._memory = '20G'

    def end(self):
        super(TableStandardAgent, self).end()


class TableStandardTool(Tool):
    def __init__(self, config):
        super(TableStandardTool, self).__init__(config)
        self.program = {
            'python': 'program/Python/bin/python'
        }
        self.script = {
            'table_kit': os.path.join(self.config.PACKAGE_DIR, 'tool_lab/table_kit.py')
        }
        self.file = {
            'transposed_table': os.path.join(self.output_dir, '{}'.format(os.path.basename(self.option("table").prop['path']))),
        }

    def run(self):
        super(TableStandardTool, self).run()

        self.run_table_kit_standard()
        self.set_output()
        self.end()

    def run_table_kit_standard(self):
        cmd = '{} {}'.format(self.program['python'], self.script['table_kit'])
        cmd += ' {}'.format('standard')
        cmd += ' --input {}'.format(self.option('table').path)
        cmd += ' --sep {}'.format(self.option('sep'))
        cmd += ' --observe {}'.format(self.option('observe'))
        if self.option('observe') == 'refrow':
            cmd += ' --guard {}'.format(self.option('guard'))
        if self.option('observe') == 'refcol':
            cmd += ' --guard {}'.format(self.option('ref_feature'))
        cmd += ' --feature {}'.format(self.option('feature'))
        cmd += ' --output {}'.format(self.file['transposed_table'])
        cmd_name = 'run_table_kit_standard'
        runcmd(self,cmd_name,cmd)

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
            'id': 'table_kit_standard_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'tool',
            'name': 'tool_lab.table_kit_standard',
            'instant': False,
            'options': {
                'table': '/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/table.txt',
                'sep': 'tab',
                'observe': 'sum',
                # 'guard': '',
                'feature': 'standard_scale',

            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)



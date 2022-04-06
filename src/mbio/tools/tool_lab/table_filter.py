# -*- coding: utf-8 -*-
# __author__ = 'zoujiaxun,qinjincheng'

import os
import unittest
import pandas as pd

from biocluster.agent import Agent
from biocluster.tool import Tool
from mbio.packages.whole_transcriptome.utils import runcmd
from biocluster.core.exceptions import OptionError


class TableFilterAgent(Agent):
    '''
    last_modify: 2019.12.12
    '''

    def __init__(self, parent):
        super(TableFilterAgent, self).__init__(parent)
        options = [
            {'name': 'table', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            # {'name': 'manipulating', 'type': 'string'},
            {'name': 'sep', 'type': 'string'},
            {'name': 'table_head', 'type': 'bool'},
            {'name': 'contain', 'type': 'bool'},
            {'name': 'lower', 'type': 'bool'},
            {'name': 'number', 'type': 'int'},
            {'name': 'info', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'filter_table', 'type': 'outfile', 'format': 'whole_transcriptome.common'}

        ]
        self.add_option(options)

    def check_options(self):
        pass

    def set_resource(self):
        self._cpu = 10
        self._memory = '20G'

    def end(self):
        super(TableFilterAgent, self).end()


class TableFilterTool(Tool):
    def __init__(self, config):
        super(TableFilterTool, self).__init__(config)
        self.program = {
            'python': 'miniconda2/bin/python'
        }
        self.script = {
            'table_kit': os.path.join(self.config.PACKAGE_DIR, 'tool_lab/table_kit.py')
        }
        self.file = {
            'filter_table': os.path.join(self.output_dir, '{}'.format(os.path.basename(self.option("table").path))),
        }

    def run(self):
        super(TableFilterTool, self).run()

        self.run_table_kit_filter()
        self.set_output()
        self.end()

    def run_table_kit_filter(self):
        cmd = '{} {}'.format(self.program['python'], self.script['table_kit'])
        cmd += ' {}'.format('filter')
        cmd += ' --input {}'.format(self.option('table').path)
        cmd += ' --sep {}'.format(self.option('sep'))
        cmd += ' --table_head {}'.format(self.option('table_head'))
        cmd += ' --contain {}'.format(self.option('contain'))
        cmd += ' --lower {}'.format(self.option('lower'))
        cmd += ' --number {}'.format(self.option('number'))
        cmd += ' --info {}'.format(self.option('info').path)
        cmd += ' --output {}'.format(self.file['filter_table'])
        cmd_name = 'run_table_kit_filter'
        runcmd(self,cmd_name,cmd)

    def set_output(self):
        self.option('filter_table').set_path(self.file['filter_table'])

class TestFunction(unittest.TestCase):
    '''
    This is test for the tool. Just run this script to do test.
    '''

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'table_kit_filter_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'tool',
            'name': 'tool_lab.table_filter',
            'instant': False,
            'options': {
                'table': '/mnt/ilustre/users/sanger-dev/workspace/20200429/TableFillna_p8qq9kcm6kr3n6h6hi6k2o6u43_0429135143916001_5840/remote_input/table/b40b752d50f52c2976aec945fea71d32.txt',
                'sep': 'tab',
                'info': '/mnt/ilustre/users/sanger-dev/workspace/20200429/TableFillna_p8qq9kcm6kr3n6h6hi6k2o6u43_0429135143916001_5840/remote_input/table/info',
                'table_head': True,
                'contain': True,
                'number':1



            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)



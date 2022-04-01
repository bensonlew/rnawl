# -*- coding: utf-8 -*-
# __author__ = 'zoujiaxun'

import os
import unittest
import pandas as pd

from biocluster.agent import Agent
from biocluster.tool import Tool
from mbio.packages.whole_transcriptome.utils import runcmd
from biocluster.core.exceptions import OptionError


class BsjreadsAgent(Agent):
    '''
    last_modify: 2019.9.17
    '''

    def __init__(self, parent):
        super(BsjreadsAgent, self).__init__(parent)
        options = [
            {'name': 'sample', 'type': 'string', 'default': None},
            {'name': 'type', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'BSJ', 'type': 'outfile', 'format': 'whole_transcriptome.common'}

        ]
        self.add_option(options)

    def check_options(self):
        pass

    def set_resource(self):
        self._cpu = 1
        self._memory = '5G'

    def end(self):
        super(BsjreadsAgent, self).end()


class BsjreadsTool(Tool):
    def __init__(self, config):
        super(BsjreadsTool, self).__init__(config)
        self.program = {
            'python': 'program/Python/bin/python'
        }
        self.script = {
            'bsjreads': os.path.join(self.config.PACKAGE_DIR, 'whole_transcriptome/circrna/bsjreads.py')
        }
        self.file = {
            'BSJ': os.path.join(self.output_dir, '{}_circRNA_merge_signal_type_bsj.txt'.format(self.option('sample'))),

        }

    def run(self):
        super(BsjreadsTool, self).run()
        self.run_bsjreads()
        self.set_output()
        self.end()

    def run_bsjreads(self):
        cmd = '{} {}'.format(self.program['python'], self.script['bsjreads'])
        cmd += ' -i {} -o {}'.format(self.option('type').path, self.file['BSJ'])
        cmd_name = 'run_bsjreads'
        runcmd(self,cmd_name,cmd)

    def set_output(self):
        self.option('BSJ').set_path(self.file['BSJ'])

class TestFunction(unittest.TestCase):
    '''
    This is test for the tool. Just run this script to do test.
    '''

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'bsjreads_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'tool',
            'name': 'whole_transcriptome.circrna.bsjreads',
            'instant': False,
            'options': {
                'type': '/mnt/ilustre/users/sanger-dev/workspace/20190924/Single_signal_1746_8723/Signal/output/circ_merge_signal.txt',
                'sample': 'zjx'

            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)



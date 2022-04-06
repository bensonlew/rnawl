# -*- coding: utf-8 -*-
# __author__ = 'zoujiaxun'

import os
import unittest
import pandas as pd

from biocluster.agent import Agent
from biocluster.tool import Tool
from mbio.packages.whole_transcriptome.utils import runcmd
from biocluster.core.exceptions import OptionError


class CirctypeAgent(Agent):
    '''
    last_modify: 2019.9.17
    '''

    def __init__(self, parent):
        super(CirctypeAgent, self).__init__(parent)
        options = [
            {'name': 'sample', 'type': 'string', 'default': None},
            {'name': 'annotate', 'type': 'infile', 'format': 'whole_transcriptome.gtf'},
            {'name': 'signal', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'type', 'type': 'outfile', 'format': 'whole_transcriptome.common'}


        ]
        self.add_option(options)

    def check_options(self):
        pass

    def set_resource(self):
        self._cpu = 1
        self._memory = '10G'

    def end(self):
        super(CirctypeAgent, self).end()


class CirctypeTool(Tool):
    def __init__(self, config):
        super(CirctypeTool, self).__init__(config)
        self.program = {
            'python': 'miniconda2/bin/python'
        }
        self.script = {
            'circtype': os.path.join(self.config.PACKAGE_DIR, 'whole_transcriptome/circrna/circRNAtype.py')
        }
        self.file = {
            'type': os.path.join(self.output_dir, '{}_circRNA_merge_signal_type.txt'.format(self.option('sample'))),

        }

    def run(self):
        super(CirctypeTool, self).run()
        self.run_type()
        self.set_output()
        self.end()

    def run_type(self):
        cmd = '{} {}'.format(self.program['python'], self.script['circtype'])
        cmd += ' -g {} -s {} -o {}'.format(self.option('annotate').path, self.option('signal').path, self.file['type'])
        cmd_name = 'run_type'
        runcmd(self,cmd_name,cmd)

    def set_output(self):
        self.option('type').set_path(self.file['type'])


class TestFunction(unittest.TestCase):
    '''
    This is test for the tool. Just run this script to do test.
    '''

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'circtype_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'tool',
            'name': 'whole_transcriptome.circrna.circtype',
            'instant': False,
            'options': {
                'annotate': '/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/data/Homo_sapiens.GRCh38.96.gtf',
                'signal': '/mnt/ilustre/users/sanger-dev/workspace/20190924/Single_signal_1746_8723/Signal/output/circ_merge_signal.txt'

            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)



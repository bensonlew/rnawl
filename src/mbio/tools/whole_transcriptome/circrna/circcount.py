# -*- coding: utf-8 -*-
# __author__ = 'zoujiaxun'

import os
import unittest
import re
import pandas as pd

from biocluster.agent import Agent
from biocluster.tool import Tool
from mbio.packages.whole_transcriptome.utils import runcmd
from biocluster.core.exceptions import OptionError


class CirccountAgent(Agent):
    '''
    last_modify: 2019.9.17
    '''

    def __init__(self, parent):
        super(CirccountAgent, self).__init__(parent)
        options = [
            {'name': 'list_id', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'counts', 'type': 'outfile', 'format': 'whole_transcriptome.common'},


        ]
        self.add_option(options)

    def check_options(self):
        pass

    def set_resource(self):
        self._cpu = 1
        self._memory = '10G'

    def end(self):
        super(CirccountAgent, self).end()


class CirccountTool(Tool):
    def __init__(self, config):
        super(CirccountTool, self).__init__(config)
        self.program = {
            'python': 'program/Python/bin/python',
        }
        self.script = {
            'circcount': os.path.join(self.config.PACKAGE_DIR, 'whole_transcriptome/circrna/circcount.py')
        }

        self.file = {
            'counts': os.path.join(self.output_dir, 'count.txt'),

        }

    def run(self):
        super(CirccountTool, self).run()
        self.run_circ_count()
        self.set_output()
        self.end()

    def run_circ_count(self):
        cmd = '{} {}'.format(self.program['python'], self.script['circcount'])
        cmd += ' -i {} -o {}'.format(self.option('list_id').path,self.file['counts'])
        cmd_name = 'run_circ_count'
        runcmd(self,cmd_name,cmd)

    def set_output(self):
        self.option('counts').set_path(self.file['counts'])


class TestFunction(unittest.TestCase):
    '''
    This is test for the tool. Just run this script to do test.
    '''

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'rpm_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'tool',
            'name': 'whole_transcriptome.circrna.circcount',
            'instant': False,
            'options': {
                'list_id': '/mnt/ilustre/users/sanger-dev/workspace/20191219/Circrna_tsg_36604/CircBrush/circ_id.list',

            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)



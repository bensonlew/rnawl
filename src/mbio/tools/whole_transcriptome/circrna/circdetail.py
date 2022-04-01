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


class CircdetailAgent(Agent):
    '''
    last_modify: 2019.9.17
    '''

    def __init__(self, parent):
        super(CircdetailAgent, self).__init__(parent)
        options = [
            {'name': 'list_id', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'details', 'type': 'outfile', 'format': 'whole_transcriptome.common'},


        ]
        self.add_option(options)

    def check_options(self):
        pass

    def set_resource(self):
        self._cpu = 1
        self._memory = '10G'

    def end(self):
        super(CircdetailAgent, self).end()


class CircdetailTool(Tool):
    def __init__(self, config):
        super(CircdetailTool, self).__init__(config)
        self.program = {
            'python': 'program/Python/bin/python',
        }
        self.script = {
            'circdetail': os.path.join(self.config.PACKAGE_DIR, 'whole_transcriptome/circrna/detail.py')
        }

        self.file = {
            'details': os.path.join(self.output_dir, 'detail.txt'),

        }

    def run(self):
        super(CircdetailTool, self).run()
        self.run_circ_detail()
        self.set_output()
        self.end()

    def run_circ_detail(self):
        cmd = '{} {}'.format(self.program['python'], self.script['circdetail'])
        cmd += ' -i {} -o {}'.format(self.option('list_id').path,self.file['details'])
        cmd_name = 'run_circ_detail'
        runcmd(self,cmd_name,cmd)

    def set_output(self):
        self.option('details').set_path(self.file['details'])


class TestFunction(unittest.TestCase):
    '''
    This is test for the tool. Just run this script to do test.
    '''

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'detail_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'tool',
            'name': 'whole_transcriptome.circrna.circdetail',
            'instant': False,
            'options': {
                'list_id': '/mnt/ilustre/users/sanger-dev/workspace/20191213/Longrna_tsg_36545/CircBrush/circ_id.list'
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)



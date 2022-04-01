# -*- coding: utf-8 -*-
# __author__ = 'zoujiaxun'

import os
import unittest
import pandas as pd

from biocluster.agent import Agent
from biocluster.tool import Tool
from mbio.packages.whole_transcriptome.utils import runcmd
from biocluster.core.exceptions import OptionError


class CutfastaAgent(Agent):
    '''
    last_modify: 2019.9.17
    '''

    def __init__(self, parent):
        super(CutfastaAgent, self).__init__(parent)
        options = [
            {'name': 'details', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'circfasta', 'type': 'infile', 'format': 'whole_transcriptome.fasta'},
            {'name': 'ref', 'type': 'outfile', 'format': 'whole_transcriptome.fasta'},
            {'name': 'new', 'type': 'outfile', 'format': 'whole_transcriptome.fasta'}

        ]
        self.add_option(options)

    def check_options(self):
        pass

    def set_resource(self):
        self._cpu = 1
        self._memory = '10G'

    def end(self):
        super(CutfastaAgent, self).end()



class CutfastaTool(Tool):
    def __init__(self, config):
        super(CutfastaTool, self).__init__(config)
        self.program = {
            'python': 'program/Python/bin/python'
        }
        self.script = {
            'cutfasta': os.path.join(self.config.PACKAGE_DIR, 'whole_transcriptome/circrna/cutfasta.py')
        }
        self.file = {
            'ref': os.path.join(self.output_dir, 'ref.fasta'),
            'new': os.path.join(self.output_dir, 'new.fasta')


        }

    def run(self):
        super(CutfastaTool, self).run()
        self.run_cutfasta()
        self.set_output()
        self.end()


    def run_cutfasta(self):
        cmd = '{} {}'.format(self.program['python'], self.script['cutfasta'])
        cmd += ' -d {} -f {} -r {} -n {}'.format(self.option('details').path, self.option('circfasta').path, self.file['ref'], self.file['new'])
        cmd_name = 'run_cutfasta'
        runcmd(self,cmd_name,cmd)

    def set_output(self):
        self.option('ref').set_path(self.file['ref'])
        self.option('new').set_path(self.file['new'])

class TestFunction(unittest.TestCase):
    '''
    This is test for the tool. Just run this script to do test.
    '''

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'cutfasta_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'tool',
            'name': 'whole_transcriptome.circrna.cutfasta',
            'instant': False,
            'options': {
                'details':'/mnt/ilustre/users/sanger-dev/workspace/20191127/Circrna_tsg_36283/CircBrush/Ifdatabase/output/detail.txt',
                'circfasta': '/mnt/ilustre/users/sanger-dev/workspace/20191127/Circrna_tsg_36283/CircBrush/Getfasta/output/circrna.fasta'

            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)



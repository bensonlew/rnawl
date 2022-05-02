# -*- coding: utf-8 -*-
# __author__ = 'zoujiaxun'
import os
import unittest
import pandas as pd

from biocluster.agent import Agent
from biocluster.tool import Tool
from mbio.packages.whole_transcriptome.utils import runcmd
from biocluster.core.exceptions import OptionError


class BlastnAgent(Agent):
    '''
    last_modify: 2019.10.21
    '''

    def __init__(self, parent):
        super(BlastnAgent, self).__init__(parent)
        options = [
            {'name': 'circfasta', 'type': 'infile', 'format': 'whole_transcriptome.fasta'},
            {'name': 'genome', 'type': 'infile', 'format': 'whole_transcriptome.fasta'},
            {'name': 'circ_blast', 'type': 'outfile', 'format': 'whole_transcriptome.common'}




        ]
        self.add_option(options)

    def check_options(self):
        pass

    def set_resource(self):
        self._cpu = 10
        self._memory = '40G'

    def end(self):
        super(BlastnAgent, self).end()


class BlastnTool(Tool):
    def __init__(self, config):
        super(BlastnTool, self).__init__(config)
        self.program = {
            # 'makeblastdb': 'miniconda2/bin/makeblastdb',
            # 'blastn': 'miniconda2/bin/blastn'
            'blat': 'bioinfo/align/ucsc_tools/blat'
        }
        self.file = {
            'circ_blat': os.path.join(self.output_dir, 'circ_blat.txt'),


        }


    def run(self):
        super(BlastnTool, self).run()
        # self.run_makeblastdb()
        self.run_blat()
        self.set_output()
        self.end()

    # def run_makeblastdb(self):
    #     cmd = '{} -in {} -dbtype nucl -parse_seqids'.format(self.program['makeblastdb'], self.option('genome_new').path)
    #     cmd_name = 'run_makeblastdb'
    #     runcmd(self, cmd_name, cmd)

    def run_blat(self):
        cmd = '{} {} {} -out=blast8 -minIdentity=100 {}'.format(self.program['blat'], self.option('genome').path, self.option('circfasta').path, self.file['circ_blat'])
        cmd_name = 'run_blastn'
        runcmd(self,cmd_name,cmd)

    def set_output(self):
        self.option('circ_blat').set_path(self.file['circ_blat'])





class TestFunction(unittest.TestCase):
    '''
    This is test for the tool. Just run this script to do test.
    '''

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'blastn_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'tool',
            'name': 'whole_transcriptome.circrna.blastn',
            'instant': False,
            'options': {
                'genome': '/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/BLASTcirc/hg19.fa',
                'circfasta': '/mnt/ilustre/users/sanger-dev/workspace/20191021/Single_getfasta_8458_3962/Getfasta/output/circrna.fasta'

            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)



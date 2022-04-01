# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import os
import unittest

from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError
from biocluster.tool import Tool

from mbio.packages.whole_transcriptome.utils import runcmd


class BwaindexAgent(Agent):
    '''
    last_modify: 2019.10.10
    '''

    def __init__(self, parent):
        super(BwaindexAgent, self).__init__(parent)
        options = [
            {'name': 'genome', 'type': 'infile', 'format': 'whole_transcriptome.fasta'},
        ]
        self.add_option(options)

    def check_options(self):
        for k, v in self._options.items():
            self.logger.debug('{} -> {}'.format(k, v.value))
        else:
            return True

    def set_resource(self):
        self._cpu = 1
        self._memory = '{}G'.format(int(os.path.getsize(self.option('genome').path) / 1024.0 ** 3 + 40))

    def end(self):
        super(BwaindexAgent, self).end()


class BwaindexTool(Tool):
    def __init__(self, config):
        super(BwaindexTool, self).__init__(config)
        self.program = {
            'bwa': 'install_packages/bwa-0.7.16a/bwa',
        }

    def run(self):
        super(BwaindexTool, self).run()
        if self.check_no_index():
            self.run_bwa_index()
        self.set_output()
        self.end()

    def check_no_index(self):
        for suffix in ('amb', 'ann', 'bwt', 'pac', 'sa'):
            if not os.path.isfile('{}.{}'.format(self.option('genome').path, suffix)):
                return True
        else:
            return False

    def run_bwa_index(self):
        cmd = '{} {} {} {} {}'.format(self.program['bwa'], 'index', '-a', 'bwtsw', self.option('genome').path)
        cmd_name = 'run_bwa_index'
        runcmd(self, cmd_name, cmd)

    def set_output(self):
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
            'id': 'bwaindex_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'tool',
            'name': 'whole_transcriptome.circrna.bwaindex',
            'instant': False,
            'options': {
                'genome': '/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/GRCh38_Ensembl_96/dna/Homo_sapiens.GRCh38.dna.toplevel.fa',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)

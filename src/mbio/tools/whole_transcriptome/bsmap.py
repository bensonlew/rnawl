# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import os
import unittest

from biocluster.agent import Agent
from biocluster.tool import Tool
from mbio.packages.whole_transcriptome.utils import runcmd


class BsmapAgent(Agent):
    """
    last_modify: 2020.02.18
    """

    def __init__(self, parent):
        super(BsmapAgent, self).__init__(parent)
        options = [
            {'name': 'query_a', 'type': 'infile', 'format': 'whole_transcriptome.fastq'},
            {'name': 'query_b', 'type': 'infile', 'format': 'whole_transcriptome.fastq'},
            {'name': 'reference', 'type': 'infile', 'format': 'whole_transcriptome.fasta'},
            {'name': 'sample', 'type': 'string', 'default': None},
            {'name': 'n_processors', 'type': 'int', 'default': 8},
            {'name': 'set_mapping_strand_information', 'type': 'bool', 'default': True},
            {'name': 'report_unmapped', 'type': 'bool', 'default': False},
            {'name': 'only_unique', 'type': 'bool', 'default': True},
            {'name': 'only_paired', 'type': 'bool', 'default': True},
            {'name': 'report_zero_meth', 'type': 'bool', 'default': True},
        ]
        self.add_option(options)

    def check_options(self):
        for k, v in self._options.items():
            self.logger.debug('{} => {}'.format(k, v.value))

    def set_resource(self):
        self._cpu = 1
        self._memory = '20G'

    def end(self):
        super(BsmapAgent, self).end()


class BsmapTool(Tool):
    def __init__(self, config):
        super(BsmapTool, self).__init__(config)
        self.program = {
            'bsmap': 'bioinfo/wgbs/miniconda3/bin/bsmap',
            'python': 'program/Python/bin/python',
        }
        self.script = {
            'methratio': os.path.join(self.config.SOFTWARE_DIR, 'bioinfo/wgbs/miniconda3/bin/methratio.py')
        }
        self.dir = {
            'result': os.path.join(self.output_dir, 'result'),
            'report': os.path.join(self.output_dir, 'report'),
        }
        self.file = {
            'bsp': os.path.join(self.output_dir, '{}.bsp'.format(self.option('sample'))),
            'ratio': os.path.join(self.output_dir, '{}.txt'.format(self.option('sample'))),
        }
        self.set_environ(PATH=os.path.join(self.config.SOFTWARE_DIR, 'bioinfo/wgbs/miniconda3/bin'))

    def run(self):
        super(BsmapTool, self).run()
        self.run_bsmap()
        self.run_methratio()
        self.set_output()
        self.end()

    def run_bsmap(self):
        cmd = '{} -a {} -b {}'.format(self.program['bsmap'], self.option('query_a').path, self.option('query_b').path)
        cmd += ' -d {}'.format(self.option('reference').path)
        cmd += ' -o {}'.format(self.file['bsp'])
        cmd += ' -p {}'.format(self.option('n_processors'))
        if self.option('set_mapping_strand_information'):
            cmd += ' -n 1'
        if self.option('report_unmapped'):
            cmd += ' -u'
        runcmd(self, 'run_bsmap', cmd)

    def run_methratio(self):
        cmd = '{} {}'.format(self.program['python'], self.script['methratio'])
        cmd += ' -o {}'.format(self.file['ratio'])
        cmd += ' -d {}'.format(self.option('reference').path)
        if self.option('only_unique'):
            cmd += ' -u'
        if self.option('only_paired'):
            cmd += ' -p'
        if self.option('report_zero_meth'):
            cmd += ' -z'
        cmd += ' {}'.format(self.file['bsp'])
        runcmd(self, 'run_methratio', cmd)

    def set_output(self):
        pass


class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'bsmap_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'tool',
            'name': 'whole_transcriptome.bsmap',
            'instant': False,
            'options': {
                'query_a': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/medip/clean_data/YGXLTHZ1225-4'
                           '-YG201902876-JJH_L2.clean.1.fastq',
                'query_b': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/medip/clean_data/YGXLTHZ1225-4'
                           '-YG201902876-JJH_L2.clean.2.fastq',
                'reference': '/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens'
                             '/Mitochondrion/rCRS/rCRS.fasta',
                'sample': 'YGXLTHZ1225-4-YG201902876-JJH_L2',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)

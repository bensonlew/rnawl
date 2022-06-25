# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import json
import os
import unittest

from biocluster.agent import Agent
from biocluster.tool import Tool

from mbio.packages.whole_transcriptome.utils import runcmd


class KallistoAgent(Agent):
    '''
    last_modify: 2019.11.19
    '''

    def __init__(self, parent):
        super(KallistoAgent, self).__init__(parent)
        options = [
            {'name': 'index', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'bias', 'type': 'bool', 'default': True},
            {'name': 'single', 'type': 'bool', 'default': False},
            {'name': 'reads', 'type': 'infile', 'format': 'whole_transcriptome.fastq'},
            {'name': 'reads_1', 'type': 'infile', 'format': 'whole_transcriptome.fastq'},
            {'name': 'reads_2', 'type': 'infile', 'format': 'whole_transcriptome.fastq'},
            {'name': 'strand_specific', 'type': 'bool', 'default': True},
            {'name': 'strand_dir', 'type': 'string', 'default': 'RF'},
            {'name': 'threads', 'type': 'int', 'default': 8},
            {'name': 't2g', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'sample', 'type': 'string', 'default': None},

            {'name': 't_abundance', 'type': 'outfile', 'format': 'whole_transcriptome.common'},
            {'name': 'g_abundance', 'type': 'outfile', 'format': 'whole_transcriptome.common'},
        ]
        self.add_option(options)

    def check_options(self):
        for k, v in self._options.items():
            self.logger.debug('{} -> {}'.format(k, v.value))
        else:
            return True

    def set_resource(self):
        self._cpu = 10
        if self.option('single'):
            byte = os.path.getsize(self.option('reads').path)
        else:
            byte = os.path.getsize(self.option('reads_1').path) + os.path.getsize(self.option('reads_2').path)
        self._memory = '{}G'.format(int(byte / 1024.0 ** 3 + 20))

    def end(self):
        super(KallistoAgent, self).end()


class KallistoTool(Tool):
    def __init__(self, config):
        super(KallistoTool, self).__init__(config)
        python_path = self.config.SOFTWARE_DIR + '/miniconda2/bin/'
        self.set_environ(PATH=python_path)
        self.program = {
            'kallisto': 'bioinfo/rna/kallisto-0.46.1/kallisto',
            'python': 'miniconda2/bin/python'
        }
        self.script = {
            'abundance': os.path.join(self.config.PACKAGE_DIR, 'whole_transcriptome/expression/abundance.py')
        }
        self.file = {
            'i_abundance': os.path.join(self.work_dir, 'abundance.tsv'),
            't_abundance': os.path.join(self.output_dir, 'T.abundance.txt'),
            'g_abundance': os.path.join(self.output_dir, 'G.abundance.txt'),
            'json': os.path.join(self.work_dir, 'input.json')
        }

    def run(self):
        super(KallistoTool, self).run()
        self.run_kallisto_quant()
        self.pre_abundance()
        self.run_abundance()
        self.set_output()
        self.end()

    def run_kallisto_quant(self):
        cmd = '{} quant'.format(self.program['kallisto'])
        cmd += ' -i {}'.format(self.option('index').path)
        cmd += ' -o {}'.format(self.work_dir)
        if self.option('bias'):
            cmd += ' --bias'
        if self.option('strand_specific'):
            if self.option('strand_dir').startswith('R'):
                cmd += ' --rf-stranded'
            elif self.option('strand_dir').startswith('F'):
                cmd += ' --fr-stranded'
        cmd += ' -t {}'.format(self.option('threads'))
        if self.option('single'):
            cmd += ' --single {}'.format(self.option('reads').path)
        else:
            cmd += ' {} {}'.format(self.option('reads_1').path, self.option('reads_2').path)
        runcmd(self, 'run_kallisto_quant', cmd)

    def pre_abundance(self):
        json.dump({
            'm_table': self.option('t2g').path,
            'i_table': self.file['i_abundance'],
            't_table': self.file['t_abundance'],
            'g_table': self.file['g_abundance']
        }, open(self.file['json'], 'w'), indent=4)

    def run_abundance(self):
        cmd = '{} {} {}'.format(self.program['python'], self.script['abundance'], self.file['json'])
        runcmd(self, 'run_abundance', cmd)

    def set_output(self):
        self.option('t_abundance').set_path(self.file['t_abundance'])
        self.option('g_abundance').set_path(self.file['g_abundance'])


class TestFunction(unittest.TestCase):
    '''
    This is test for the tool. Just run this script to do test.
    '''

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'kallisto_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'tool',
            'name': 'whole_transcriptome.expression.kallisto',
            'instant': False,
            'options': {
                'index': '/mnt/ilustre/users/sanger-dev/workspace/20191119/Single_index_1614_1657/Index/output/kallisto.index',
                'reads_1': '/mnt/ilustre/users/sanger-dev/workspace/20191105/Longrna_tsg_36051/FastpRna/fastq_dir/C_6w.clean.1.fastq',
                'reads_2': '/mnt/ilustre/users/sanger-dev/workspace/20191105/Longrna_tsg_36051/FastpRna/fastq_dir/C_6w.clean.2.fastq',
                't2g': '/mnt/ilustre/users/sanger-dev/workspace/20191105/Longrna_tsg_36051/LargeGush/MergeKnownNew/t2g.txt',
                'sample': 'C_6w',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)

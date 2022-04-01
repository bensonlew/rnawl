# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import os
import unittest

from biocluster.agent import Agent
from biocluster.tool import Tool

from mbio.packages.whole_transcriptome.utils import runcmd


class IndexAgent(Agent):
    '''
    last_modify: 2019.11.19
    '''

    def __init__(self, parent):
        super(IndexAgent, self).__init__(parent)
        PROGRAM = ('rsem', 'kallisto', 'salmon')
        options = [
            {'name': 'program', 'type': 'string', 'default': PROGRAM[2]},
            {'name': 'transcripts', 'type': 'infile', 'format': 'whole_transcriptome.fasta'},
            {'name': 't2g', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'kmer_size', 'type': 'int', 'default': 31},
            {'name': 'kmer_len', 'type': 'int', 'default': 31},
            {'name': 'threads', 'type': 'int', 'default': 8},
            {'name': 'rsem_index', 'type': 'outfile', 'format': 'whole_transcriptome.common'},
            {'name': 'kallisto_index', 'type': 'outfile', 'format': 'whole_transcriptome.common'},
            {'name': 'salmon_index', 'type': 'outfile', 'format': 'whole_transcriptome.common_dir'},
            {'name': 'perfect_hash', 'type': 'bool', 'default': True},
        ]
        self.add_option(options)

    def check_options(self):
        for k, v in self._options.items():
            self.logger.debug('{} -> {}'.format(k, v.value))
        else:
            return True

    def set_resource(self):
        self._cpu = 1
        self._memory = '5G'

    def end(self):
        super(IndexAgent, self).end()


class IndexTool(Tool):
    def __init__(self, config):
        super(IndexTool, self).__init__(config)
        python_path = self.config.SOFTWARE_DIR + '/program/Python/bin/'
        self.set_environ(PATH=python_path)
        self.program = {
            'bowtie2': 'bioinfo/rna/miniconda3/bin/bowtie2',
            'rsem_prepare_reference': 'bioinfo/rna/RSEM-1.3.1/rsem-prepare-reference',
            'kallisto': 'bioinfo/rna/kallisto-0.46.1/kallisto',
            'salmon': 'bioinfo/rna/salmon-0.14.1/bin/salmon',
        }
        self.file = {
            'g2t': os.path.join(self.work_dir, 'knownIsoforms.txt'),
            'index': os.path.join(self.output_dir, '{}.index'.format(self.option('program')))
        }
        self.dir = {
            'bowtie2_path': os.path.join(self.config.SOFTWARE_DIR, 'bioinfo/rna/miniconda3/bin'),
            'perl_path': os.path.join(self.config.SOFTWARE_DIR, 'program/perl/perls/perl-5.24.0/bin')
        }
        self.set_environ(PATH=self.dir['perl_path'])

    def run(self):
        super(IndexTool, self).run()
        if self.option('program') == 'rsem':
            self.run_rsem_prepare_referece()
        elif self.option('program') == 'kallisto':
            self.run_kallisto_index()
        elif self.option('program') == 'salmon':
            self.run_salmon_index()
        self.set_output()
        self.end()

    def run_rsem_prepare_referece(self):
        open(self.file['g2t'], 'w').writelines(
            '{}\t{}\n'.format(*line.strip().split('\t')[::-1]) for line in open(self.option('t2g').path))
        cmd = '{}'.format(self.program['rsem_prepare_reference'])
        cmd += ' --transcript-to-gene-map {}'.format(self.file['g2t'])
        cmd += ' --bowtie2 --bowtie2-path {}'.format(self.dir['bowtie2_path'])
        cmd += ' --num-threads {}'.format(self.option('threads'))
        cmd += ' {} {}'.format(self.option('transcripts').path, self.file['index'])
        runcmd(self, 'run_rsem_prepare_referece', cmd)

    def run_kallisto_index(self):
        cmd = '{} index'.format(self.program['kallisto'])
        cmd += ' -i {}'.format(self.file['index'])
        cmd += ' -k {}'.format(self.option('kmer_size'))
        cmd += ' {}'.format(self.option('transcripts').path)
        runcmd(self, 'run_kallisto_index', cmd)

    def run_salmon_index(self):
        cmd = '{} index'.format(self.program['salmon'])
        cmd += ' -t {}'.format(self.option('transcripts').path)
        cmd += ' -k {}'.format(self.option('kmer_len'))
        cmd += ' -i {}'.format(self.file['index'])
        cmd += ' -p {}'.format(self.option('threads'))
        if self.option('perfect_hash'):
            cmd += ' --perfectHash'
        runcmd(self, 'run_salmon_index', cmd)

    def set_output(self):
        if self.option('program') == 'rsem':
            open(self.file['index'], 'w').close()
            self.option('rsem_index').set_path(self.file['index'])
        elif self.option('program') == 'kallisto':
            self.option('kallisto_index').set_path(self.file['index'])
        elif self.option('program') == 'salmon':
            self.option('salmon_index').set_path(self.file['index'])


class TestFunction(unittest.TestCase):
    '''
    This is test for the tool. Just run this script to do test.
    '''

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'index_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'tool',
            'name': 'whole_transcriptome.expression.index',
            'instant': False,
            'options': {
                'program': 'kallisto',
                'transcripts': '/mnt/ilustre/users/sanger-dev/workspace/20191105/Longrna_tsg_36051/LargeGush/MergeKnownNew/output/known_and_new.fa',
                't2g': '/mnt/ilustre/users/sanger-dev/workspace/20191105/Longrna_tsg_36051/LargeGush/MergeKnownNew/output/t2g.txt',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)

# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import os
import unittest

from biocluster.agent import Agent
from biocluster.tool import Tool
from mbio.packages.whole_transcriptome.utils import runcmd


class Hisat2Agent(Agent):
    """
    last_modify: 2020.03.04
    """

    def __init__(self, parent):
        super(Hisat2Agent, self).__init__(parent)
        options = [
            {'name': 'ht2_idx', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'is_se', 'type': 'bool', 'default': False},
            {'name': 'm1', 'type': 'infile', 'format': 'whole_transcriptome.fastq'},
            {'name': 'm2', 'type': 'infile', 'format': 'whole_transcriptome.fastq'},
            {'name': 'r', 'type': 'infile', 'format': 'whole_transcriptome.fastq'},
            {'name': 'sample_name', 'type': 'string', 'default': None},
            {'name': 'threads', 'type': 'int', 'default': 8},
            {'name': 'output_fmt', 'type': 'string', 'default': 'BAM'},
            {'name': 'bam', 'type': 'outfile', 'format': 'whole_transcriptome.common'},
        ]
        self.add_option(options)

    def check_options(self):
        for k, v in self._options.items():
            self.logger.debug('{} => {}'.format(k, v.value))

    def set_resource(self):
        self._cpu = 1
        self._memory = '20G'

    def end(self):
        super(Hisat2Agent, self).end()


class Hisat2Tool(Tool):
    def __init__(self, config):
        super(Hisat2Tool, self).__init__(config)
        self.program = {
            'hisat2': 'bioinfo/align/hisat2/hisat2-2.1.0/hisat2',
            'samtools': 'bioinfo/rna/miniconda2/bin/samtools'
        }
        self.file = {
            'sam': os.path.join(self.work_dir, '{}.sam'.format(self.option('sample_name'))),
            'bam': os.path.join(self.output_dir, '{}.bam'.format(self.option('sample_name')))
        }

    def run(self):
        super(Hisat2Tool, self).run()
        self.run_hisat2()
        self.run_samtools()
        self.set_output()
        self.end()

    def run_hisat2(self):
        cmd = '{} -x {}'.format(self.program['hisat2'], self.option('ht2_idx').path)
        if self.option('is_se'):
            cmd += ' -U {}'.format(self.option('r').path)
        else:
            cmd += ' -1 {}'.format(self.option('m1').path)
            cmd += ' -2 {}'.format(self.option('m2').path)
        cmd += ' -S {}'.format(self.file['sam'])
        cmd += ' -p {}'.format(self.option('threads'))
        runcmd(self, 'run_hisat2', cmd)

    def run_samtools(self):
        cmd = '{} sort -o {}'.format(self.program['samtools'], self.file['bam'])
        cmd += ' -O {}'.format(self.option('output_fmt'))
        cmd += ' -@ {}'.format(self.option('threads'))
        cmd += ' {}'.format(self.file['sam'])
        runcmd(self, 'run_samtools', cmd)

    def set_output(self):
        self.option('bam').set_path(self.file['bam'])


class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'hisat2_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'tool',
            'name': 'whole_transcriptome.hisat2',
            'instant': False,
            'options': {
                'ht2_idx': '/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/bacilli'
                           '/Geobacillus_thermodenitrificans/GCF_000015745.1_ASM1574v1/dna/GCF_000015745'
                           '.1_ASM1574v1_genomic',
                'm1': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/data/Geobacillus_thermodenitrificans'
                      '/clean_data/zt_005.clean.1.fastq',
                'm2': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/data/Geobacillus_thermodenitrificans'
                      '/clean_data/zt_005.clean.2.fastq',
                'sample_name': 'zt_005',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)

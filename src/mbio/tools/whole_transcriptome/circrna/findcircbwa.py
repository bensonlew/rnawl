# -*- coding: utf-8 -*-
# __author__ = 'zoujiaxun'

import os
import unittest

import pandas as pd
from biocluster.agent import Agent
from biocluster.tool import Tool

from mbio.packages.whole_transcriptome.utils import runcmd


class FindcircbwaAgent(Agent):
    '''
    last_modify: 2019.8.22
    '''

    def __init__(self, parent):
        super(FindcircbwaAgent, self).__init__(parent)
        options = [
            {'name': 'sample', 'type': 'string', 'default': None},
            {'name': 'genome', 'type': 'infile', 'format': 'whole_transcriptome.fasta'},
            {'name': 'find_circ', 'type': 'outfile', 'format': 'whole_transcriptome.common'},
            {'name': 'rnasam', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'rnabam', 'type': 'outfile', 'format': 'whole_transcriptome.common'},
            {'name': 'find_circ_new', 'type': 'outfile', 'format': 'whole_transcriptome.common'},
            {'name': 'find_circ_circ_splice_sites', 'type': 'outfile', 'format': 'whole_transcriptome.common'},
            {'name': 'find_circ_filter', 'type': 'outfile', 'format': 'whole_transcriptome.common'},
            {'name': 'junction_reads', 'type': 'int', 'default': 2},
            {'name': 'circrna_length', 'type': 'int', 'default': 100000}
        ]
        self.add_option(options)
        self._memory_increase_step = 30

    def check_options(self):
        pass

    def set_resource(self):
        self._cpu = 1
        self._memory = '30G'

    def end(self):
        super(FindcircbwaAgent, self).end()


class FindcircbwaTool(Tool):
    def __init__(self, config):
        super(FindcircbwaTool, self).__init__(config)

        self.program = {
            'python': 'bioinfo/rna/miniconda2/bin/python',
            'samtools': 'bioinfo/align/samtools-1.3.1/samtools',
        }
        self.file = {
            'find_circ': os.path.join(self.output_dir, '{}_find_circ2'.format(self.option('sample'))),
            'rnabam': os.path.join(self.output_dir, '{}.bam'.format(self.option('sample'))),
            'find_circ_circ_splice_sites': os.path.join(self.output_dir, '{}_find_circ2/circ_splice_sites.bed'.format(
                self.option('sample'))),
            'find_circ_new': os.path.join(self.output_dir,
                                          '{}_find_circ2/{}_circ_splice_sites.bed'.format(self.option('sample'),
                                                                                          self.option('sample'))),
            'find_circ_filter': os.path.join(self.output_dir, '{}_find_circ_filter.txt'.format(self.option('sample')))

        }
        self.script = {
            'find_circ2': os.path.join(self.config.SOFTWARE_DIR, 'bioinfo/circ_rna/find_circ2-master/find_circ.py'),
        }

    def run(self):
        super(FindcircbwaTool, self).run()
        self.run_bam()
        self.run_findcirc2()
        self.run_rename()
        self.filter()
        self.set_output()
        self.end()

    def run_bam(self):
        cmd = '{} view -hbuS {} -o {}'.format(self.program['samtools'], self.option('rnasam').path, self.file['rnabam'])
        cmd_name = 'run_bam'
        runcmd(self, cmd_name, cmd)

    def run_findcirc2(self):
        cmd = '{} {}'.format(self.program['python'], self.script['find_circ2'])
        cmd += ' -G {}'.format(self.option('genome').path)
        cmd += ' -a 20'
        cmd += ' -o {}'.format(self.file['find_circ'])
        cmd += ' {}'.format(self.file['rnabam'])
        cmd_name = 'run_findcirc2'
        runcmd(self, cmd_name, cmd)

    def run_rename(self):
        os.rename(self.file['find_circ_circ_splice_sites'], self.file['find_circ_new'])
        # cmd = 'mv {} {}'.format(self.file['find_circ_circ_splice_sites'], self.file['find_circ_new'])
        # cmd.software_dir = "/usr/bin/"
        # cmd_name = 'run_rename'
        # runcmd(self, cmd_name, cmd, shell=True)

    def filter(self):
        f1 = pd.read_table(self.file['find_circ_new'])
        f2 = f1[f1['n_frags'] >= self.option('junction_reads')]
        f3 = f2[~f2['#chrom'].isin(['MT'])]
        f4 = f3[f3['category'].isnull()]
        f5 = f4[f4['end'] - f4['start'] < self.option('circrna_length')]
        f6 = f5.rename(
            columns={'#chrom': 'chr', 'start': 'circRNA_start', 'end': 'circRNA_end', 'n_frags': 'junction_reads',
                     'name': 'circRNA_name'})
        f7 = f6[['circRNA_name', 'chr', 'strand', 'circRNA_start', 'circRNA_end', 'junction_reads']]
        f7['circRNA_start'] = f7['circRNA_start'] + 1
        index = f7.index.tolist()
        f7['circRNA_name'] = ['{}_{}_{}'.format(f7['chr'][i], f7['circRNA_start'][i], f7['circRNA_end'][i]) for i in
                              index]
        f7.sort_values(by=['chr', 'circRNA_start', 'circRNA_end'], ascending=True, inplace=True)
        f7.to_csv(self.file['find_circ_filter'], index=False, sep='\t')

    def set_output(self):
        self.option('find_circ_filter').set_path(self.file['find_circ_filter'])


class TestFunction(unittest.TestCase):
    '''
    This is test for the tool. Just run this script to do test.
    '''

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'findcircbwa_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'tool',
            'name': 'whole_transcriptome.circrna.findcircbwa',
            'instant': False,
            'options': {
                'genome': '/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/GRCh38_Ensembl_96/dna/Homo_sapiens.GRCh38.dna.toplevel.fa',
                'rnasam': '/mnt/ilustre/users/sanger-dev/workspace/20191014/Longrna_workflow_9038_1820/CircBrush/Ciriaddfindcirc4/Bwamem/output/C3.sam',
                'sample': 'C3'
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)

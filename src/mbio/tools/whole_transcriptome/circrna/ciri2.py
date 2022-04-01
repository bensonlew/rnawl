# -*- coding: utf-8 -*-
# __author__ = 'zoujiaxun'

import os
import unittest

import pandas as pd
from biocluster.agent import Agent
from biocluster.tool import Tool

from mbio.packages.whole_transcriptome.utils import runcmd


class Ciri2Agent(Agent):
    '''
    last_modify: 2019.8.22
    '''

    def __init__(self, parent):
        super(Ciri2Agent, self).__init__(parent)
        options = [
            {'name': 'sample', 'type': 'string', 'default': None},
            {'name': 'genome', 'type': 'infile', 'format': 'whole_transcriptome.fasta'},
            # {'name': 'fq1', 'type': 'infile', 'format': 'ref_rna_v2.fastq'},
            # {'name': 'fq2', 'type': 'infile', 'format': 'ref_rna_v2.fastq'},
            {'name': 'annotate', 'type': 'infile', 'format': 'whole_transcriptome.gtf'},
            {'name': 'rnasam', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'circrna', 'type': 'outfile', 'format': 'whole_transcriptome.common'},
            {'name': 'ciri_circ_filter', 'type': 'outfile', 'format': 'whole_transcriptome.common'},
            {'name': 'max_span', 'type': 'int', 'default': 200000},
            {'name': 'strigrncy', 'type': 'string', 'default': 'high'},
            {'name': 'mapq_uni', 'type': 'int', 'default': 10},
            {'name': 'chrM', 'type': 'string', 'default': 'chrM'},
            {'name': 'num_thread', 'type': 'int', 'default': 1},
            {'name': 'junction_reads', 'type': 'int', 'default': 2},
            {'name': 'circrna_length', 'type': 'int', 'default': 100000}
            # {'name': 'rnasam', 'type': 'infile', 'format': 'ref_rna_v2.common'}
        ]
        self.add_option(options)
        self._memory_increase_step = 100

    def check_options(self):
        pass

    def set_resource(self):
        self._cpu = 2
        self._memory = '20G'

    def end(self):
        super(Ciri2Agent, self).end()


class Ciri2Tool(Tool):
    def __init__(self, config):
        super(Ciri2Tool, self).__init__(config)

        self.program = {
            'bwa': 'install_packages/bwa-0.7.16a/bwa',
            'python': 'program/Python/bin/python',
            'perl': 'program/perl-5.24.0/bin/perl'
        }
        self.file = {
            'circrna': os.path.join(self.output_dir, '{}.txt'.format(self.option('sample'))),
            'ciri_circ_filter': os.path.join(self.output_dir, '{}_ciri_circ_filter.txt'.format(self.option('sample')))
        }
        self.script = {
            'ciri2': os.path.join(self.config.SOFTWARE_DIR, 'bioinfo/circ_rna/CIRI_v2.0.6/CIRI_v2.0.6/CIRI2.pl')
        }

    def run(self):
        super(Ciri2Tool, self).run()
        self.run_ciri2()
        self.filter()
        self.set_output()
        self.end()

    def run_ciri2(self):
        cmd = '{}'.format(self.program['perl'])
        cmd += ' {}'.format(self.script['ciri2'])
        cmd += ' -I {}'.format(self.option('rnasam').path)

        cmd += ' -O {}'.format(self.file['circrna'])
        cmd += ' -F {}'.format(self.option('genome').path)
        cmd += ' -S {}'.format(self.option('max_span'))
        cmd += ' -{}'.format(self.option('strigrncy'))
        cmd += ' -U {}'.format(self.option('mapq_uni'))
        cmd += ' -M {}'.format(self.option('chrM'))
        cmd += ' -T {}'.format(self.option('num_thread'))
        cmd += ' -A {}'.format(self.option('annotate').path)
        cmd_name = 'run_ciri2'
        command = self.add_command(cmd_name, cmd, ignore_error=True)
        command.run()
        self.wait()
        if command.return_code == 0:
            pass
        elif command.return_code in [-9,-11]:
            self.logger.info("return code: {}".format(command.return_code))
            self.add_state('memory_limit', 'memory is low!')
        elif command.return_code is None:
            command.rerun()
            if command.return_code == 0:
                self.logger.info("succeed in running run_ciri2")
                pass
        else:
            self.set_error("fail to run run_ciri2")
        # runcmd(self, cmd_name, cmd)

    def filter(self):
        f1 = pd.read_table(self.file['circrna'])
        f2 = f1[f1['circRNA_end'] - f1['circRNA_start'] < self.option('circrna_length')]
        f3 = f2[f2['#junction_reads'] >= self.option('junction_reads')]
        f4 = f3.rename(columns={'circRNA_ID': 'circRNA_name', '#junction_reads': 'junction_reads'})
        f5 = f4[['circRNA_name', 'chr', 'strand', 'circRNA_start', 'circRNA_end', 'junction_reads']]
        index = f5.index.tolist()
        f5['circRNA_name'] = ['{}_{}_{}'.format(f5['chr'][i], f5['circRNA_start'][i], f5['circRNA_end'][i]) for i in
                              index]
        f5.to_csv(self.file['ciri_circ_filter'], index=False, sep='\t')

    def set_output(self):
        self.option('ciri_circ_filter').set_path(self.file['ciri_circ_filter'])


class TestFunction(unittest.TestCase):
    '''
    This is test for the tool. Just run this script to do test.
    '''

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'ciri2_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'tool',
            'name': 'whole_transcriptome.circrna.ciri2',
            'instant': False,
            'options': {
                'genome': '/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/data/Homo_sapiens.GRCh38.dna.toplevel.fa',
                'max_span': 200000,
                'strigrncy': 'high',
                'mapq_uni': 10,
                'chrM': 'chrM',
                'num_thread': 1,
                'annotate': '/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/data/Homo_sapiens.GRCh38.96.gtf',
                'rnasam': '/mnt/ilustre/users/sanger-dev/workspace/20190929/Single_circrna_1983_3415/Circrna/Bwamem__1/output/zjx.sam',
                'sample': 'zjx'

            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)

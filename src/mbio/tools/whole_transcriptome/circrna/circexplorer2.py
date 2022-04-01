# -*- coding: utf-8 -*-
# __author__ = 'zoujiaxun'

import os
import unittest
import pandas as pd

from biocluster.agent import Agent
from biocluster.tool import Tool
from mbio.packages.whole_transcriptome.utils import runcmd
from biocluster.core.exceptions import OptionError


class Circexplorer2Agent(Agent):
    '''
    last_modify: 2019.9.29
    '''

    def __init__(self, parent):
        super(Circexplorer2Agent, self).__init__(parent)
        options = [
            {'name': 'genome', 'type': 'infile', 'format': 'whole_transcriptome.fasta'},
            # {'name': 'annotate', 'type': 'infile', 'format': 'whole_transcriptome.gtf'},
            # {'name': 'fq1', 'type': 'infile', 'format': 'whole_transcriptome.fastq'},
            # {'name': 'fq2', 'type': 'infile', 'format': 'whole_transcriptome.fastq'},
            {'name': 'GenePred_ref', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'rnasam', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name':'circexplorer2_circRNA', 'type': 'outfile', 'format': 'whole_transcriptome.common'},
            {'name': 'sample', 'type': 'string', 'default': None},
            {'name': 'circexplorer_filter', 'type': 'outfile', 'format': 'whole_transcriptome.common'},
            {'name': 'junction_reads', 'type': 'int', 'default': 2},
            {'name': 'circrna_length', 'type': 'int', 'default': 100000}

        ]
        self.add_option(options)

    def check_options(self):
        pass

    def set_resource(self):
        self._cpu = 10
        self._memory = '50G'

    def end(self):
        super(Circexplorer2Agent, self).end()


class Circexplorer2Tool(Tool):
    def __init__(self, config):
        super(Circexplorer2Tool, self).__init__(config)
        self.program = {
            'Circexplorer2': 'bioinfo/rna/miniconda2/bin/CIRCexplorer2',
            'bwa': 'install_packages/bwa-0.7.16a/bwa',
            'python': 'program/Python/bin/python',
            'perl': 'program/perl-5.24.0/bin/perl',


        }
        self.file = {
            'circexplorer2_circRNA': os.path.join(self.output_dir, '{}_circexplorer2_circRNA.txt'.format(self.option('sample'))),
            # 'rnasam': os.path.join(self.output_dir, '{}_rnasam.sam'.format(self.option('sample'))),
            'circexplorer_filter': os.path.join(self.output_dir, '{}_circexplorer_filter.txt'.format(self.option('sample')))




        }


    def run(self):
        super(Circexplorer2Tool, self).run()
        # self.run_bwa_index()
        # self.run_bwa_mem()
        self.run_circexplore_parse()
        self.run_circexplore_annotate()
        self.run_filter()
        self.set_output()
        self.end()




    # def run_bwa_index(self):
    #     cmd = '{} {} {} {} {}'.format(self.program['bwa'], 'index', '-a', 'bwtsw', self.option('genome').path)
    #     cmd_name = 'run_bwa_index'
    #     runcmd(self, cmd_name, cmd)
    #
    # def run_bwa_mem(self):
    #     cmd = '{} {} {} {} {} {} {} {} {}'.format(self.program['bwa'], 'mem', '-T', 19, self.option('genome').path, self.option('fq1').path, self.option('fq2').path, '-o', self.file['rnasam'])
    #     cmd_name = 'run_bwa_mem'
    #     runcmd(self, cmd_name, cmd)

    def run_circexplore_parse(self):
        cmd = '{} parse'.format(self.program['Circexplorer2'])
        cmd += ' -t BWA'
        cmd += ' {}'.format(self.option('rnasam').path)
        cmd += ' >Circexplore_parse.log'
        cmd_name = 'run_circexplore_parse'
        runcmd(self, cmd_name, cmd, shell=True)

    def run_circexplore_annotate(self):
        cmd = '{} annotate'.format(self.program['Circexplorer2'])
        cmd += ' -r {}'.format(self.option('GenePred_ref').path)
        cmd += ' -g {}'.format(self.option('genome').path)
        cmd += ' -b back_spliced_junction.bed'
        cmd += ' -o {}'.format(self.file['circexplorer2_circRNA'])
        cmd += ' > CIRIexplorer2_annotate.log'
        cmd_name = 'run_circexplore_annotate'
        runcmd(self, cmd_name, cmd, shell=True)

    def run_filter(self):
        f1 = pd.read_table(self.file['circexplorer2_circRNA'],header = None)
        f2 = f1[[0,1,2,5,12]]
        f3 = f2.rename(columns={0: 'chr', 1: 'circRNA_start', 2: 'circRNA_end', 5: 'strand', 12: 'junction_reads'})
        f4 = f3[f3['circRNA_end'] - f3['circRNA_start'] < self.option('circrna_length')]
        f5 = f4[f4['junction_reads'] >= self.option('junction_reads')]
        f5['circRNA_start'] = f5['circRNA_start'] + 1
        index = f5.index.tolist()
        f5['circRNA_name'] = ['{}_{}_{}'.format(f5['chr'][i], f5['circRNA_start'][i], f5['circRNA_end'][i]) for i in index]
        f5.sort_values(by=['chr', 'circRNA_start', 'circRNA_end'], ascending=True, inplace=True)
        f5.to_csv(self.file['circexplorer_filter'],index = False,sep = '\t')


    def set_output(self):
        self.option('circexplorer_filter').set_path(self.file['circexplorer_filter'])


class TestFunction(unittest.TestCase):
    '''
    This is test for the tool. Just run this script to do test.
    '''

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'circexplorer2_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'tool',
            'name': 'whole_transcriptome.circrna.circexplorer2',
            'instant': False,
            'options': {

                'genome': '/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/data/Homo_sapiens.GRCh38.dna.toplevel.fa',
                'GenePred_ref': '/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/data/Gh38_ref.txt',
                # 'fq1': '/mnt/ilustre/users/isanger/workspace/20190213/Single_LncrnaQc_7789/LncrnaQc/output/sickle_dir/Con1_sickle_l.fastq',
                # 'fq2': '/mnt/ilustre/users/isanger/workspace/20190213/Single_LncrnaQc_7789/LncrnaQc/output/sickle_dir/Con1_sickle_r.fastq',
                'sample': 'zjx',
                'rnasam': '/mnt/ilustre/users/sanger-dev/workspace/20190929/Single_circrna_1983_3415/Circrna/Bwamem__1/output/zjx.sam'


            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)



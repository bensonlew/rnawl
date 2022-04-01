# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import os
import shutil
import unittest

from biocluster.agent import Agent
from biocluster.tool import Tool

from mbio.packages.denovo_rna_v2.functions import reminder, runcmd


class DenovoUploadAgent(Agent):
    '''
    last_modify: 2019.08.05
    '''

    def __init__(self, parent):
        PROGRAM_OPTIONS = ('blastn', 'blastx', 'tblastx')
        DBTYPE_OPTIONS = ('nucl', 'prot')
        EVALUE_OPTIONS = (10, 1e-1)
        WORD_SIZE_OPTIONS = (11, 3)
        super(DenovoUploadAgent, self).__init__(parent)
        options = [
            {'name': 'query', 'type': 'infile', 'format': 'denovo_rna_v2.fasta'},
            {'name': 'subject', 'type': 'infile', 'format': 'denovo_rna_v2.fasta'},
            {'name': 'program', 'type': 'string', 'default': PROGRAM_OPTIONS[0]},
            {'name': 'dbtype', 'type': 'string', 'default': DBTYPE_OPTIONS[0]},
            {'name': 'evalue', 'type': 'float', 'default': EVALUE_OPTIONS[0]},
            {'name': 'word_size', 'type': 'int', 'default': WORD_SIZE_OPTIONS[0]},
            {'name': 'max_hsps', 'type': 'int', 'default': 100},
            {'name': 'outfmt', 'type': 'string', 'default': '6'},
            {'name': 'num_threads', 'type': 'int', 'default': 8},
            {'name': 'tabular', 'type': 'outfile', 'format': 'denovo_rna_v2.common'},
        ]
        self.add_option(options)

    @reminder
    def check_options(self):
        for k, v in self._options.items():
            self.logger.debug('{} = {}'.format(k, v.value))

    @reminder
    def set_resource(self):
        self._cpu = 1
        self._memory = '{}G'.format(int(os.path.getsize(self.option('subject').path) / 1024.0 ** 3 * 2 + 20))

    @reminder
    def end(self):
        super(DenovoUploadAgent, self).end()


class DenovoUploadTool(Tool):
    def __init__(self, config):
        super(DenovoUploadTool, self).__init__(config)
        self.program = {
            'makeblastdb': 'bioinfo/denovo_rna_v2/miniconda2/bin/makeblastdb',
            'blastn': 'bioinfo/denovo_rna_v2/miniconda2/bin/blastn',
            'blastp': 'bioinfo/denovo_rna_v2/miniconda2/bin/blastp',
            'blastx': 'bioinfo/denovo_rna_v2/miniconda2/bin/blastx',
            'tblastn': 'bioinfo/denovo_rna_v2/miniconda2/bin/tblastn',
            'tblastx': 'bioinfo/denovo_rna_v2/miniconda2/bin/tblastx'
        }
        self.file = {
            'database': os.path.join(self.work_dir, 'subject.fasta'),
            'out': os.path.join(self.output_dir, 'blast.m8.txt')
        }

    @reminder
    def run(self):
        super(DenovoUploadTool, self).run()
        self.run_makeblastdb()
        {'blastn': self.run_blastn, 'blastp': self.run_blastp, 'blastx': self.run_blastx, 'tblastn': self.run_tblastn,
         'tblastx': self.run_tblastx}[self.option('program')]()
        self.set_output()
        self.end()

    @reminder
    def run_makeblastdb(self):
        shutil.copy(self.option('subject').path, self.file['database'])
        cmd = '{} -in {} -dbtype {}'.format(self.program['makeblastdb'], self.file['database'], self.option('dbtype'))
        cmd_name = 'run_makeblastdb'
        runcmd(self, cmd_name, cmd)

    @reminder
    def run_blastn(self):
        cmd = '{} -db {}'.format(self.program['blastn'], self.file['database'])
        cmd += ' -query {}'.format(self.option('query').path)
        cmd += ' -out {}'.format(self.file['out'])
        cmd += ' -evalue {}'.format(self.option('evalue'))
        cmd += ' -word_size {}'.format(self.option('word_size'))
        cmd += ' -max_hsps {}'.format(self.option('max_hsps'))
        cmd += ' -outfmt {}'.format(self.option('outfmt'))
        cmd += ' -num_threads {}'.format(self.option('num_threads'))
        cmd_name = 'run_blastn'
        runcmd(self, cmd_name, cmd)

    @reminder
    def run_blastp(self):
        cmd = '{} -db {}'.format(self.program['blastp'], self.file['database'])
        cmd += ' -query {}'.format(self.option('query').path)
        cmd += ' -out {}'.format(self.file['out'])
        cmd += ' -evalue {}'.format(self.option('evalue'))
        cmd += ' -word_size {}'.format(self.option('word_size'))
        cmd += ' -max_hsps {}'.format(self.option('max_hsps'))
        cmd += ' -outfmt {}'.format(self.option('outfmt'))
        cmd += ' -num_threads {}'.format(self.option('num_threads'))
        cmd_name = 'run_blastp'
        runcmd(self, cmd_name, cmd)

    @reminder
    def run_blastx(self):
        cmd = '{} -db {}'.format(self.program['blastx'], self.file['database'])
        cmd += ' -query {}'.format(self.option('query').path)
        cmd += ' -out {}'.format(self.file['out'])
        cmd += ' -evalue {}'.format(self.option('evalue'))
        cmd += ' -word_size {}'.format(self.option('word_size'))
        cmd += ' -max_hsps {}'.format(self.option('max_hsps'))
        cmd += ' -outfmt {}'.format(self.option('outfmt'))
        cmd += ' -num_threads {}'.format(self.option('num_threads'))
        cmd_name = 'run_blastn'
        runcmd(self, cmd_name, cmd)

    @reminder
    def run_tblastn(self):
        cmd = '{} -db {}'.format(self.program['tblastn'], self.file['database'])
        cmd += ' -query {}'.format(self.option('query').path)
        cmd += ' -out {}'.format(self.file['out'])
        cmd += ' -evalue {}'.format(self.option('evalue'))
        cmd += ' -word_size {}'.format(self.option('word_size'))
        cmd += ' -max_hsps {}'.format(self.option('max_hsps'))
        cmd += ' -outfmt {}'.format(self.option('outfmt'))
        cmd += ' -num_threads {}'.format(self.option('num_threads'))
        cmd_name = 'run_tblastxn'
        runcmd(self, cmd_name, cmd)

    @reminder
    def run_tblastx(self):
        cmd = '{} -db {}'.format(self.program['tblastx'], self.file['database'])
        cmd += ' -query {}'.format(self.option('query').path)
        cmd += ' -out {}'.format(self.file['out'])
        cmd += ' -evalue {}'.format(self.option('evalue'))
        cmd += ' -word_size {}'.format(self.option('word_size'))
        cmd += ' -max_hsps {}'.format(self.option('max_hsps'))
        cmd += ' -outfmt {}'.format(self.option('outfmt'))
        cmd += ' -num_threads {}'.format(self.option('num_threads'))
        cmd_name = 'run_tblastx'
        runcmd(self, cmd_name, cmd)

    @reminder
    def set_output(self):
        self.option('tabular').set_path(self.file['out'])


class TestFunction(unittest.TestCase):
    '''
    This is test for the tool. Just run this script to do test.
    '''

    def test_blastn(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'denovo_upload_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'tool',
            'name': 'denovo_upl_v2.alignment.denovo_upload',
            'instant': False,
            'options': {
                'query': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/denovo_rna_v2/Trinity.fasta',
                'subject': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/denovo_rna_v2/transcript.fa',
                'program': 'blastn',
                'evalue': 10,
                'word_size': 11,
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

    def test_tblastx(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'denovo_upload_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'tool',
            'name': 'denovo_rna_v2.alignment.denovo_upload',
            'instant': False,
            'options': {
                'query': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/denovo_rna_v2/Trinity.fasta',
                'subject': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/denovo_rna_v2/cds.fa',
                'program': 'tblastx',
                'evalue': 1e-1,
                'word_size': 3,
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

    def test_blastx(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'denovo_upload_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'tool',
            'name': 'denovo_rna_v2.alignment.denovo_upload',
            'instant': False,
            'options': {
                'query': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/denovo_rna_v2/Trinity.fasta',
                'subject': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/denovo_rna_v2/peptide.fa',
                'program': 'blastx',
                'dbtype': 'prot',
                'evalue': 1e-1,
                'word_size': 3,
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test_blastn'),
                    TestFunction('test_tblastx'),
                    TestFunction('test_blastx')])
    unittest.TextTestRunner(verbosity=2).run(suite)

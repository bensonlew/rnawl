# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import os
import shutil
import unittest

from biocluster.agent import Agent
from biocluster.config import Config
from biocluster.tool import Tool

from mbio.packages.denovo_rna_v2.functions import reminder, runcmd


class AlignFinalAgent(Agent):
    '''
    last_modify: 2019.08.09
    '''

    def __init__(self, parent):
        super(AlignFinalAgent, self).__init__(parent)
        ABOUT_OPTIONS = ('genome', 'database', 'upload')
        DB_TYPE_OPTIONS = ('transcript', 'cds', 'peptide')
        options = [
            {'name': 'tabular', 'type': 'infile', 'format': 'denovo_rna_v2.common'},
            {'name': 'genome_id', 'type': 'string', 'default': None},
            {'name': 'evalue', 'type': 'float', 'default': 1e-1},
            {'name': 'about', 'type': 'string', 'default': ABOUT_OPTIONS[0]},
            {'name': 'db_type', 'type': 'string', 'default': ABOUT_OPTIONS[0]},
            {'name': 'result', 'type': 'outfile', 'format': 'denovo_rna_v2.common'},
        ]
        self.add_option(options)

    @reminder
    def check_options(self):
        for k, v in self._options.items():
            self.logger.debug('{} = {}'.format(k, v.value))

    @reminder
    def set_resource(self):
        self._cpu = 1
        self._memory = '8G'

    @reminder
    def end(self):
        super(AlignFinalAgent, self).end()


class AlignFinalTool(Tool):
    def __init__(self, config):
        super(AlignFinalTool, self).__init__(config)
        self.program = {
            'python': 'program/Python/bin/python',
            'gtftogenepred': 'bioinfo/align/ucsc_tools/gtfToGenePred',
            'genepredtobed': 'bioinfo/align/ucsc_tools/genePredToBed',
            'bedtools': 'bioinfo/denovo_rna_v2/miniconda2/bin/bedtools',
        }
        self.script = {
            'get_bed6': os.path.join(self.config.PACKAGE_DIR, 'denovo_rna_v2/alignment/get_bed6.py'),
            'get_result': os.path.join(self.config.PACKAGE_DIR, 'denovo_rna_v2/alignment/get_result.py'),
        }
        self.file = {
            'result': os.path.join(self.output_dir, 'result.{}.txt'.format(self.option('about')))
        }
        self.genome_doc = Config().get_mongo_client(mtype='ref_rna_v2', dydb_forbid=True)[Config().get_mongo_dbname('ref_rna_v2', dydb_forbid=True)][
            'sg_genome_db'].find_one({'genome_id': self.option('genome_id')})

    @reminder
    def run(self):
        super(AlignFinalTool, self).run()
        if self.option('about') == 'genome':
            self.run_gtftogenepred()
            self.run_genepredtobed()
            self.run_get_bed6()
            self.run_bedtools()
        self.run_get_result()
        self.set_output()
        self.end()

    @reminder
    def run_gtftogenepred(self):
        self.file['gtf'] = os.path.join(self.work_dir, '{}.gtf'.format(self.option('genome_id')))
        self.file['genepred'] = os.path.join(self.work_dir, '{}.genepred'.format(self.option('genome_id')))
        shutil.copy(os.path.join(self.config.SOFTWARE_DIR, 'database/Genome_DB_finish', self.genome_doc['gtf']),
                    self.file['gtf'])
        cmd = '{} {} {}'.format(self.program['gtftogenepred'], self.file['gtf'], self.file['genepred'])
        cmd_name = 'run_gtftogenepred'
        runcmd(self, cmd_name, cmd)

    @reminder
    def run_genepredtobed(self):
        self.file['ref_bed'] = os.path.join(self.work_dir, '{}.bed'.format(self.option('genome_id')))
        cmd = '{} {} {}'.format(self.program['genepredtobed'], self.file['genepred'], self.file['ref_bed'])
        cmd_name = 'run_genepredtobed'
        runcmd(self, cmd_name, cmd)

    @reminder
    def run_get_bed6(self):
        self.file['hit_bed'] = os.path.join(self.work_dir, 'hits.6.bed')
        cmd = '{} {}'.format(self.program['python'], self.script['get_bed6'])
        cmd += ' --input {}'.format(self.option('tabular').path)
        cmd += ' --evalue {}'.format(self.option('evalue'))
        cmd += ' --output {}'.format(self.file['hit_bed'])
        cmd_name = 'run_get_bed6'
        runcmd(self, cmd_name, cmd)

    @reminder
    def run_bedtools(self):
        self.file['intersect'] = os.path.join(self.work_dir, 'intersect.txt')
        cmd = '{} intersect -s -wo'.format(self.program['bedtools'])
        cmd += ' -a {}'.format(self.file['hit_bed'])
        cmd += ' -b {}'.format(self.file['ref_bed'])
        cmd += ' > {}'.format(self.file['intersect'])
        cmd_name = 'run_bedtools'
        runcmd(self, cmd_name, cmd, shell=True)

    def run_get_result(self):
        if self.option('about') != 'upload':
            self.file['biomart'] = os.path.join(self.work_dir, '{}.biomart.txt'.format(self.option('genome_id')))
            self.file['g2t2p'] = os.path.join(self.work_dir, '{}.g2t2p.txt'.format(self.option('genome_id')))
            shutil.copy(os.path.join(self.config.SOFTWARE_DIR, 'database/Genome_DB_finish', self.genome_doc['bio_mart_annot']),
                        self.file['biomart'])
            shutil.copy(os.path.join(self.config.SOFTWARE_DIR, 'database/Genome_DB_finish', self.genome_doc['g2t2p']),
                        self.file['g2t2p'])
        cmd = '{} {} {}'.format(self.program['python'], self.script['get_result'], self.option('about'))
        cmd += ' --m8_out {}'.format(self.option('tabular').path)
        if self.option('about') == 'genome':
            cmd += ' --bt_out {}'.format(self.file['intersect'])
            cmd += ' --bm_file {}'.format(self.file['biomart'])
            cmd += ' --bm_type {}'.format(self.genome_doc['biomart_gene_annotype'])
        elif self.option('about') == 'database':
            cmd += ' --db_type {}'.format(self.option('db_type'))
            cmd += ' --bm_file {}'.format(self.file['biomart'])
            cmd += ' --bm_type {}'.format(self.genome_doc['biomart_gene_annotype'])
            cmd += ' --g2t2p {}'.format(self.file['g2t2p'])
        cmd += ' --output {}'.format(self.file['result'])
        cmd_name = 'run_get_result'
        runcmd(self, cmd_name, cmd)

    @reminder
    def set_output(self):
        self.option('result').set_path(self.file['result'])


class TestFunction(unittest.TestCase):
    '''
    This is test for the tool. Just run this script to do test.
    '''

    def test_genome(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'align_final_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'tool',
            'name': 'denovo_rna_v2.alignment.align_final',
            'instant': False,
            'options': {
                'tabular': '/mnt/ilustre/users/sanger-dev/workspace/20190808/Alignment_alignment_5349_7825/DenovoRef/output/blast.m8.txt',
                'genome_id': 'GM0311',
                'evalue': 1e-1,
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

    def test_database(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'align_final_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'tool',
            'name': 'denovo_rna_v2.alignment.align_final',
            'instant': False,
            'options': {
                'tabular': '/mnt/ilustre/users/sanger-dev/workspace/20190809/Alignment_alignment_8342_3145/DenovoRef/output/blast.m8.txt',
                'genome_id': 'GM0265',
                'about': 'database',
                'db_type': 'peptide',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

    def test_upload(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'align_final_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'tool',
            'name': 'denovo_rna_v2.alignment.align_final',
            'instant': False,
            'options': {
                'tabular': '/mnt/ilustre/users/sanger-dev/workspace/20190812/Alignment_alignment_6503_5186/DenovoUpload/output/blast.m8.txt',
                'about': 'upload',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test_upload')])
    unittest.TextTestRunner(verbosity=2).run(suite)

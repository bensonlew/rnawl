# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import os
import unittest

from biocluster.workflow import Workflow


class AlignmentWorkflow(Workflow):
    '''
    last_modify: 2019.08.08
    '''

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(AlignmentWorkflow, self).__init__(wsheet_object)
        QUERY_SOURCE_OPTIONS = ('denovo', 'upload')
        EXP_LEVEL_OPTIONS = ('G', 'T')
        QUERY_TYPE_OPTIONS = ('nucl', 'prot')
        SUBJECT_SOURCE_OPTIONS = ('ref', 'upload', 'denovo')
        SUBJECT_TARGET_OPTIONS = ('genome', 'transcript', 'cds', 'peptide', 'raw', 'filter')
        SUBJECT_TYPE_OPTIONS = ('nucl', 'prot')
        PROGRAM_OPTIONS = ('blat', 'blastn', 'blastx', 'blastp', 'tblastn')
        options = [
            {'name': 'query_source', 'type': 'string', 'default': QUERY_SOURCE_OPTIONS[0]},
            {'name': 'query_file', 'type': 'infile', 'format': 'denovo_rna_v2.fasta'},
            {'name': 'exp_level', 'type': 'string', 'default': EXP_LEVEL_OPTIONS[0]},
            {'name': 'geneset_id', 'type': 'string', 'default': 'All'},
            {'name': 'query_type', 'type': 'string', 'default': QUERY_TYPE_OPTIONS[0]},
            {'name': 'subject_source', 'type': 'string', 'default': SUBJECT_SOURCE_OPTIONS[0]},
            {'name': 'subject_target', 'type': 'string', 'default': SUBJECT_TARGET_OPTIONS[0]},
            {'name': 'subject_file', 'type': 'infile', 'format': 'denovo_rna_v2.fasta'},
            {'name': 'subject_type', 'type': 'string', 'default': SUBJECT_TYPE_OPTIONS[0]},
            {'name': 'genome_id', 'type': 'string', 'default': None},
            {'name': 'program', 'type': 'string', 'default': PROGRAM_OPTIONS[0]},
            {'name': 'evalue', 'type': 'float', 'default': 1e-1},
            {'name': 'max_hsps', 'type': 'int', 'default': 10},
            {'name': 'result', 'type': 'outfile', 'format': 'denovo_rna_v2.common'},
            {'name': 'main_id', 'type': 'string', 'default': None},
            {'name': 'update_info', 'type': 'string', 'default': None}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())

    def check_options(self):
        for k, v in self.sheet.options().items():
            self.logger.debug('{} = {}'.format(k, v))

    def run(self):
        self.run_align_begin()
        super(AlignmentWorkflow, self).run()

    def run_align_begin(self):
        self.align_begin = self.add_tool('denovo_rna_v2.alignment.align_begin')
        if self.option('query_source') == 'denovo':
            options = {
                'purpose': 'query',
                'unit_type': self.option('query_type'),
                'task_id': self.sheet.task_id,
                'exp_level': self.option('exp_level'),
                'geneset_id': self.option('geneset_id')
            }
        elif self.option('subject_source') == 'denovo':
            options = {
                'purpose': 'subject',
                'task_id': self.sheet.task_id,
                'source': self.option('subject_target')
            }
        self.align_begin.set_options(options)
        if self.option('query_source') == 'denovo':
            if self.option('subject_source') == 'ref':
                self.align_begin.on('end', self.run_denovo_ref)
            elif self.option('subject_source') == 'upload':
                self.align_begin.on('end', self.run_denovo_upload)
        elif self.option('query_source') == 'upload':
            self.align_begin.on('end', self.run_upload_denovo)
        self.align_begin.run()

    def run_denovo_ref(self):
        self.denovo_ref = self.add_tool('denovo_rna_v2.alignment.denovo_ref')
        self.denovo_ref.set_options({
            'query': self.align_begin.option('denovo'),
            'genome_id': self.option('genome_id'),
            'target': self.option('subject_target'),
            'program': self.option('program'),
            'dbtype': self.option('subject_type'),
            'evalue': self.option('evalue'),
            'max_hsps': self.option('max_hsps'),
            'word_size': 3 if 'prot' in (self.option('query_type'), self.option('subject_type')) else 11
        })
        self.denovo_ref.on('end', self.run_align_final, 'denovo_ref')
        self.denovo_ref.run()

    def run_denovo_upload(self):
        self.denovo_upload = self.add_tool('denovo_rna_v2.alignment.denovo_upload')
        self.denovo_upload.set_options({
            'query': self.align_begin.option('denovo'),
            'subject': self.option('subject_file'),
            'program': self.option('program'),
            'dbtype': self.option('subject_type'),
            'evalue': self.option('evalue'),
            'max_hsps': self.option('max_hsps'),
            'word_size': 3 if 'prot' in (self.option('query_type'), self.option('subject_type')) else 11
        })
        self.denovo_upload.on('end', self.run_align_final, 'denovo_upload')
        self.denovo_upload.run()

    def run_upload_denovo(self):
        self.upload_denovo = self.add_tool('denovo_rna_v2.alignment.upload_denovo')
        self.upload_denovo.set_options({
            'query': self.option('query_file'),
            'subject': self.align_begin.option('denovo'),
            'program': self.option('program'),
            'dbtype': self.option('subject_type'),
            'evalue': self.option('evalue'),
            'max_hsps': self.option('max_hsps'),
            'word_size': 3 if 'prot' in (self.option('query_type'), self.option('subject_type')) else 11
        })
        self.upload_denovo.on('end', self.run_align_final, 'upload_denovo')
        self.upload_denovo.run()

    def run_align_final(self, event):
        self.align_final = self.add_tool('denovo_rna_v2.alignment.align_final')
        if event['data'] == 'denovo_ref':
            tabular = self.denovo_ref.option('tabular').path
            if self.option('subject_target') == 'genome':
                about = 'genome'
            elif self.option('subject_target') in ('transcript', 'cds', 'peptide'):
                about = 'database'
                db_type = self.option('subject_target')
        elif event['data'] == 'denovo_upload':
            tabular = self.denovo_upload.option('tabular').path
            about = 'upload'
        elif event['data'] == 'upload_denovo':
            tabular = self.upload_denovo.option('tabular').path
            about = 'upload'
        options = {
            'tabular': tabular,
            'genome_id': self.option('genome_id'),
            'evalue': self.option('evalue'),
            'about': about,
        }
        if self.option('subject_target') in ('transcript', 'cds', 'peptide'):
            options.update({'db_type': db_type})
        self.align_final.set_options(options)
        self.align_final.on('end', self.set_output, about)
        self.align_final.run()

    def set_output(self, event):
        source = self.align_final.option('result').path
        link_name = os.path.join(self.output_dir, os.path.basename(source))
        os.link(source, link_name)
        self.option('result').set_path(link_name)
        self.set_db(event['data'])

    def set_db(self, about):
        api = self.api.api('denovo_rna_v2.alignment')
        api.add_blast_detail(self.option('result').path, about, self.option('main_id'))
        self.end()

    def end(self):
        super(AlignmentWorkflow, self).end()


class TestFunction(unittest.TestCase):
    '''
    This is test for the workflow. Just run this script to do test.
    '''

    def test_denovo_ref(self):
        from mbio.workflows.denovo_rna_v2.report.alignment import AlignmentWorkflow
        from biocluster.wsheet import Sheet
        import random
        data = {
            'id': 'alignment_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'workflow',
            'name': 'denovo_rna_v2.report.alignment',
            'options': {
                'query_source': 'denovo',
                'exp_level': 'G',
                'geneset_id': '5d37b5f5c6598d2f578b456c',
                # 'exp_level': 'T',
                # 'geneset_id': '5cb64c2a17b2bf6187b2918b',
                'query_type': 'prot',
                'subject_source': 'ref',
                'subject_target': 'transcript',
                'subject_type': 'nucl',
                'genome_id': 'GM0265',
                'program': 'tblastn',
                'evalue': 1e-1,
                'max_hsps': 10
            }
        }
        wsheet = Sheet(data=data)
        wf = AlignmentWorkflow(wsheet)
        wf.sheet.id = 'tsg_33857'
        wf.sheet.project_sn = '188_5bc7e35556cba'
        wf.sheet.task_id = 'tsg_33857'
        wf.IMPORT_REPORT_DATA = False
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.run()

    def test_denovo_upload(self):
        from mbio.workflows.denovo_rna_v2.report.alignment import AlignmentWorkflow
        from biocluster.wsheet import Sheet
        import random
        data = {
            'id': 'alignment_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'workflow',
            'name': 'denovo_rna_v2.report.alignment',
            'options': {
                'query_source': 'denovo',
                'exp_level': 'T',
                'geneset_id': '5cb64c2a17b2bf6187b2918b',
                'query_type': 'prot',
                'subject_source': 'upload',
                # 'subject_file': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/denovo_rna_v2/transcript.fa',
                # 'subject_file': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/denovo_rna_v2/cds.fa',
                'subject_file': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/denovo_rna_v2/peptide.fa',
                'subject_type': 'prot',
                'program': 'blastp',
                'evalue': 1e-1,
                'max_hsps': 10
            }
        }
        wsheet = Sheet(data=data)
        wf = AlignmentWorkflow(wsheet)
        wf.sheet.id = 'tsg_33857'
        wf.sheet.project_sn = '188_5bc7e35556cba'
        wf.sheet.task_id = 'tsg_33857'
        wf.IMPORT_REPORT_DATA = False
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.run()

    def test_upload_denovo(self):
        from mbio.workflows.denovo_rna_v2.report.alignment import AlignmentWorkflow
        from biocluster.wsheet import Sheet
        import random
        data = {
            'id': 'alignment_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'workflow',
            'name': 'denovo_rna_v2.report.alignment',
            'options': {
                'query_source': 'upload',
                'query_type': 'prot',
                'query_file': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/denovo_rna_v2/peptide.fa',
                'subject_source': 'denovo',
                'subject_target': 'filter',
                'program': 'tblastn',
                'evalue': 1e-1,
                'max_hsps': 10
            }
        }
        wsheet = Sheet(data=data)
        wf = AlignmentWorkflow(wsheet)
        wf.sheet.id = 'tsg_33857'
        wf.sheet.project_sn = '188_5bc7e35556cba'
        wf.sheet.task_id = 'tsg_33857'
        wf.IMPORT_REPORT_DATA = False
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.run()

if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test_upload_denovo')])
    unittest.TextTestRunner(verbosity=2).run(suite)

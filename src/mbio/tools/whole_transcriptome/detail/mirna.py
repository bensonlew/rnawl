# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import json
import os
import unittest

from biocluster.agent import Agent
from biocluster.config import Config
from biocluster.tool import Tool
from mbio.packages.whole_transcriptome.utils import runcmd


class MirnaAgent(Agent):
    '''
    last_modify: 2019.10.31
    '''

    def __init__(self, parent):
        super(MirnaAgent, self).__init__(parent)
        options = [
            {'name': 'known_mirna_detail', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'novel_mirna_detail', 'type': 'infile', 'format': 'whole_transcriptome.common'},
        ]
        self.add_option(options)

    def check_options(self):
        for k, v in self._options.items():
            self.logger.debug('{} = {}'.format(k, v.value))

    def set_resource(self):
        self._cpu = 1
        self._memory = '5G'

    def end(self):
        super(MirnaAgent, self).end()


class MirnaTool(Tool):
    def __init__(self, config):
        super(MirnaTool, self).__init__(config)
        self.program = {
            'gffread': 'bioinfo/rna/cufflinks-2.2.1/gffread',
            'gtftogenepred': 'bioinfo/align/ucsc_tools/gtfToGenePred',
            'genepredtobed': 'bioinfo/align/ucsc_tools/genePredToBed',
            'python': 'program/Python/bin/python',
            'bedtools': 'bioinfo/ref_rna_v2/miniconda2/bin/bedtools'
        }
        self.script = {
            'get_gene_bed': os.path.join(self.config.PACKAGE_DIR, 'ref_rna_v2/get_gene_bed.py'),
            'fasta_clean': os.path.join(self.config.PACKAGE_DIR, 'whole_transcriptome/detail/fasta_clean.py'),
            'base_build': os.path.join(self.config.PACKAGE_DIR, 'whole_transcriptome/detail/mi_base_build.py')
        }
        self.file = {
            'input_json': os.path.join(self.work_dir, 'input.json'),
            'database': os.path.join(self.output_dir, 'mirna.db'),
            'gene_detail': os.path.join(self.output_dir, 'gene_detail.pk'),
            'transcript_detail': os.path.join(self.output_dir, 'transcript_detail.pk')
        }

    def run(self):
        super(MirnaTool, self).run()
        self.pre_base_build()
        self.run_base_build()
        self.set_output()
        self.end()


    def pre_base_build(self):
        json.dump({'known_mirna_detail': self.option('known_mirna_detail').path,
                   'novel_mirna_detail': self.option('novel_mirna_detail').path,
                   'transcript_detail_pk': self.file['transcript_detail']},
                  open(self.file['input_json'], 'w'), indent=4)

    def run_base_build(self):
        cmd = '{} {}'.format(self.program['python'], self.script['base_build'])
        cmd += ' --json {}'.format(self.file['input_json'])
        cmd += ' --database {}'.format(self.file['database'])
        runcmd(self, 'run_base_build', cmd)

    def set_output(self):
        pass



class TestFunction(unittest.TestCase):
    '''
    This is test for the tool. Just run this script to do test.
    '''

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'mirna_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'tool',
            'name': 'whole_transcriptome.detail.mirna',
            'instant': False,
            'options': {
                'known_mirna_detail': '/mnt/ilustre/users/sanger-dev/workspace/20191106/WholeTranscriptome_tsg_36088/Transfer1/output/srna/known_mirna/known_mirna_detail.xls',
                'novel_mirna_detail': '/mnt/ilustre/users/sanger-dev/workspace/20191106/WholeTranscriptome_tsg_36088/Transfer1/output/srna/novel_mirna/novel_mirna_detail.xls',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)

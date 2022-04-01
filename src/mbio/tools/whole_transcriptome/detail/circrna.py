# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import json
import os
import unittest

from biocluster.agent import Agent
from biocluster.config import Config
from biocluster.tool import Tool
from mbio.packages.whole_transcriptome.utils import runcmd


class CircrnaAgent(Agent):
    '''
    last_modify: 2019.10.31
    '''

    def __init__(self, parent):
        super(CircrnaAgent, self).__init__(parent)
        options = [
            {'name': 'circrna_detail', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'biomart_file', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'biomart_type', 'type': 'string', 'default': None},
        ]
        self.add_option(options)

    def check_options(self):
        for k, v in self._options.items():
            self.logger.debug('{} = {}'.format(k, v.value))

    def set_resource(self):
        self._cpu = 1
        self._memory = '5G'

    def end(self):
        super(CircrnaAgent, self).end()


class CircrnaTool(Tool):
    def __init__(self, config):
        super(CircrnaTool, self).__init__(config)
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
            'base_build': os.path.join(self.config.PACKAGE_DIR, 'whole_transcriptome/detail/circ_base_build.py')
        }
        self.file = {
            'input_json': os.path.join(self.work_dir, 'input.json'),
            'database': os.path.join(self.output_dir, 'circrna.db'),
            'gene_detail': os.path.join(self.output_dir, 'gene_detail.pk'),
            'transcript_detail': os.path.join(self.output_dir, 'transcript_detail.pk'),
            'circrna_fasta': os.path.join(os.path.dirname(self.option('circrna_detail').path), "circrna.fasta")
        }

    def run(self):
        super(CircrnaTool, self).run()
        self.pre_base_build()
        self.run_base_build()
        self.set_output()
        self.end()


    def pre_base_build(self):
        json.dump({'circrna_detail': self.option('circrna_detail').path,
                   'biomart_file': self.option('biomart_file').path,
                   'biomart_type': self.option('biomart_type'),
                   'transcript_detail_pk': self.file['transcript_detail'],
                   'circrna_fasta': self.file['circrna_fasta']},
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
            'id': 'circrna_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'tool',
            'name': 'whole_transcriptome.detail.circrna',
            'instant': False,
            'options': {
                'circrna_detail': '/mnt/ilustre/users/sanger-dev/workspace/20191106/WholeTranscriptome_tsg_36088/Transfer/output/circ_brush/detail.txt',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)

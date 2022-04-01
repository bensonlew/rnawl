# -*- coding: utf-8 -*-
# __author__ = 'zoujiaxun'

import os
import re
import unittest

import pandas as pd
from biocluster.agent import Agent
from biocluster.tool import Tool

from mbio.packages.whole_transcriptome.utils import runcmd


class GenepredAgent(Agent):
    '''
    last_modify: 2019.9.17
    '''

    def __init__(self, parent):
        super(GenepredAgent, self).__init__(parent)
        options = [
            {'name': 'annotate', 'type': 'infile', 'format': 'whole_transcriptome.gtf'},
            {'name': 'GenePred', 'type': 'outfile', 'format': 'whole_transcriptome.common'},
            {'name': 'GenePred_ref', 'type': 'outfile', 'format': 'whole_transcriptome.common'}

        ]
        self.add_option(options)

    def check_options(self):
        pass

    def set_resource(self):
        self._cpu = 10
        self._memory = '10G'

    def end(self):
        super(GenepredAgent, self).end()


class GenepredTool(Tool):
    def __init__(self, config):
        super(GenepredTool, self).__init__(config)
        self.program = {
            'gtftoGenePred': 'bioinfo/align/ucsc_tools/gtfToGenePred'
        }
        self.script = {
            'merge': os.path.join(self.config.PACKAGE_DIR, 'whole_transcriptome/circrna/merge.py')
        }
        self.file = {
            'GenePred': os.path.join(self.output_dir, 'GenePred.txt'),
            'GenePred_ref': os.path.join(self.output_dir, 'GenePred_ref.txt')

        }

    def run(self):
        super(GenepredTool, self).run()
        self.run_gtftogenepred()
        self.run_genepred_ref()
        self.set_output()
        self.end()

    def run_gtftogenepred(self):
        cmd = '{} {} {}'.format(self.program['gtftoGenePred'], self.option('annotate').path, self.file['GenePred'])
        cmd_name = 'run_gtftogenepred'
        runcmd(self, cmd_name, cmd)

    def run_genepred_ref(self):
        trans_id = []
        genes_id = []
        with open(self.option('annotate').path) as f:
            for lines in f:
                if not lines.startswith('#'):
                    (chrom, source, type, start, end, score, strand, phase, attribute_string) = lines.rstrip().split(
                        '\t')
                    tran_id = re.findall(r"transcript_id\s\"(.*?)\"", attribute_string)
                    gene_id = re.findall(r"gene_id\s\"(.*?)\"", attribute_string)
                    if tran_id:
                        trans_id.extend(tran_id)
                        genes_id.extend(gene_id)
        file = pd.DataFrame({1: genes_id, 0: trans_id})
        file_qc = file.drop_duplicates()
        gtf_txt = pd.read_table(self.file['GenePred'], header=None)
        file_merge = pd.merge(file_qc, gtf_txt, how='outer', on=0)
        file_merge = file_merge[['1_x', 0, '1_y', 2, 3, 4, 5, 6, 7, 8, 9]]
        file_merge.to_csv(self.file['GenePred_ref'], index=False, header=None, sep='\t')

    def set_output(self):
        self.option('GenePred_ref').set_path(self.file['GenePred_ref'])


class TestFunction(unittest.TestCase):
    '''
    This is test for the tool. Just run this script to do test.
    '''

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'genepred_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'tool',
            'name': 'whole_transcriptome.circrna.genepred',
            'instant': False,
            'options': {
                'annotate': '/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/data/Homo_sapiens.GRCh38.96.gtf',

            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)

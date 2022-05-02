# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import os
import unittest

from biocluster.agent import Agent
from biocluster.tool import Tool

from mbio.packages.whole_transcriptome.utils import runcmd


class StringtieAgent(Agent):
    '''
    last_modify: 2019.09.05
    '''

    def __init__(self, parent):
        super(StringtieAgent, self).__init__(parent)
        options = [
            {'name': 'input_bam', 'type': 'infile', 'format': 'align.bwa.bam'},
            {'name': 'strand_specific', 'type': 'bool', 'default': True},
            {'name': 'strand_dir', 'type': 'string', 'default': 'RF'},
            {'name': 'guide_gff', 'type': 'infile', 'format': 'gene_structure.gtf'},
            {'name': 'out_gtf', 'type': 'outfile', 'format': 'gene_structure.gtf'},
            {'name': 'cpus', 'type': 'int', 'default': 8},
            {'name': 'min_anchor_cov', 'type': 'int', 'default': 3},
            {'name': 'min_bundle_cov', 'type': 'int', 'default': 5},
            {'name': 'verbose', 'type': 'bool', 'default': True},
        ]
        self.add_option(options)

    def check_options(self):
        for k, v in self._options.items():
            self.logger.debug('{} -> {}'.format(k, v.value))
        else:
            return True

    def set_resource(self):
        self._cpu = 1
        self._memory = '{}G'.format(int(os.path.getsize(self.option('input_bam').path) / 1024.0 ** 3 * 2 + 20))

    def end(self):
        super(StringtieAgent, self).end()


class StringtieTool(Tool):
    def __init__(self, config):
        super(StringtieTool, self).__init__(config)
        self.program = {
            'stringtie': 'miniconda2/bin/stringtie'
        }
        self.file = {
            'out_gtf': os.path.join(self.output_dir,
                                    '{}.gtf'.format(os.path.basename(self.option('input_bam').path)[:-4]))
        }

    def run(self):
        super(StringtieTool, self).run()
        self.run_stringtie()
        self.set_output()
        self.end()

    def run_stringtie(self):
        cmd = '{} {}'.format(self.program['stringtie'], self.option('input_bam').path)
        if self.option('strand_specific'):
            if self.option('strand_dir').startswith('R'):
                cmd += ' --rf'
            elif self.option('strand_dir').startswith('F'):
                cmd += ' --fr'
        cmd += ' -G {}'.format(self.option('guide_gff').path)
        cmd += ' -o {}'.format(self.file['out_gtf'])
        cmd += ' -p {}'.format(self.option('cpus'))
        cmd += ' -j {}'.format(self.option('min_anchor_cov'))
        cmd += ' -c {}'.format(self.option('min_bundle_cov'))
        if self.option('verbose'):
            cmd += ' -v'
        runcmd(self, 'run_stringtie', cmd)

    def set_output(self):
        self.option('out_gtf').set_path(self.file['out_gtf'])


class TestFunction(unittest.TestCase):
    '''
    This is test for the tool. Just run this script to do test.
    '''

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'stringtie_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'tool',
            'name': 'whole_transcriptome.assembly.stringtie',
            'instant': False,
            'options': {
                'input_bam': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/data/test_bam/Con1.bam',
                'guide_gff': '/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/GRCh38_Ensembl_96/gtf/Homo_sapiens.GRCh38.96.gtf',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)

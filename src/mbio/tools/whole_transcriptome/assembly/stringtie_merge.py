# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import os
import unittest

from biocluster.agent import Agent
from biocluster.tool import Tool

from mbio.packages.whole_transcriptome.utils import runcmd


class StringtieMergeAgent(Agent):
    '''
    last_modify: 2019.09.05
    '''

    def __init__(self, parent):
        super(StringtieMergeAgent, self).__init__(parent)
        options = [
            {'name': 'gtf_list', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'guide_gff', 'type': 'infile', 'format': 'gene_structure.gtf'},
            {'name': 'out_gtf', 'type': 'outfile', 'format': 'gene_structure.gtf'},
            {'name': 'min_cov', 'type': 'int', 'default': 5},
            {'name': 'min_tpm', 'type': 'int', 'default': 1},
            {'name': 'min_iso', 'type': 'float', 'default': 0.1},
        ]
        self.add_option(options)

    def check_options(self):
        for k, v in self._options.items():
            self.logger.debug('{} -> {}'.format(k, v.value))
        else:
            return True

    def set_resource(self):
        self._cpu = 1
        self._memory = '{}G'.format(int(os.path.getsize(self.option('guide_gff').path) / 1024.0 ** 3 * 2 + 20))

    def end(self):
        super(StringtieMergeAgent, self).end()


class StringtieMergeTool(Tool):
    def __init__(self, config):
        super(StringtieMergeTool, self).__init__(config)
        self.program = {
            'stringtie': 'miniconda2/bin/stringtie'
        }
        self.file = {
            'out_gtf': os.path.join(self.output_dir, 'merged.gtf')
        }

    def run(self):
        super(StringtieMergeTool, self).run()
        self.run_stringtie_merge()
        self.set_output()
        self.end()

    def run_stringtie_merge(self):
        cmd = '{} --merge {}'.format(self.program['stringtie'], self.option('gtf_list').path)
        cmd += ' -G {}'.format(self.option('guide_gff').path)
        cmd += ' -o {}'.format(self.file['out_gtf'])
        cmd += ' -c {}'.format(self.option('min_cov'))
        cmd += ' -T {}'.format(self.option('min_tpm'))
        cmd += ' -f {}'.format(self.option('min_iso'))
        runcmd(self, 'run_stringtie_merge', cmd)

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
            'id': 'stringtie_merge_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'tool',
            'name': 'whole_transcriptome.assembly.stringtie_merge',
            'instant': False,
            'options': {
                'gtf_list': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/lnc_rna/assemble/assembly_GTF_list.txt',
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

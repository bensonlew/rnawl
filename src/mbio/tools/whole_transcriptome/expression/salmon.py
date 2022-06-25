# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import os
import unittest

from biocluster.agent import Agent
from biocluster.tool import Tool

from mbio.packages.whole_transcriptome.utils import runcmd


class SalmonAgent(Agent):
    '''
    last_modify: 2019.09.06
    '''

    def __init__(self, parent):
        super(SalmonAgent, self).__init__(parent)
        options = [
            {'name': 'lib_type', 'type': 'string', 'default': 'A'},
            {'name': 'index', 'type': 'infile', 'format': 'whole_transcriptome.common_dir'},
            {'name': 'unmated', 'type': 'infile', 'format': 'whole_transcriptome.fastq'},
            {'name': 'mates1', 'type': 'infile', 'format': 'whole_transcriptome.fastq'},
            {'name': 'mates2', 'type': 'infile', 'format': 'whole_transcriptome.fastq'},
            {'name': 'gc_bias', 'type': 'bool', 'default': True},
            {'name': 'threads', 'type': 'int', 'default': 8},
            {'name': 'gene_map', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'validate_mappings', 'type': 'bool', 'default': True},
            {'name': 'write_unmapped_names', 'type': 'bool', 'default': True},
            {'name': 'sample', 'type': 'string', 'default': None},
            {'name': 'quant_sf', 'type': 'outfile', 'format': 'whole_transcriptome.common'},
            {'name': 'quant_genes_sf', 'type': 'outfile', 'format': 'whole_transcriptome.common'},
        ]
        self.add_option(options)

    def check_options(self):
        for k, v in self._options.items():
            self.logger.debug('{} -> {}'.format(k, v.value))
        else:
            return True

    def set_resource(self):
        self._cpu = 10
        if self.option('unmated').is_set:
            memory = int(os.path.getsize(self.option('unmated').path) / 1024.0 ** 3 + 20)
        else:
            memory = int((os.path.getsize(self.option('mates1').path) + os.path.getsize(
                self.option('mates2').path)) / 1024.0 ** 3 + 20)
        self._memory = '{}G'.format(memory)

    def end(self):
        super(SalmonAgent, self).end()


class SalmonTool(Tool):
    def __init__(self, config):
        super(SalmonTool, self).__init__(config)
        python_path = self.config.SOFTWARE_DIR + '/miniconda2/bin/'
        self.set_environ(PATH=python_path)
        self.program = {
            'salmon': 'bioinfo/rna/salmon-0.14.1/bin/salmon'
        }
        self.file = {
            'quant_sf': os.path.join(self.output_dir, 'quant.sf'),
            'quant_genes_sf': os.path.join(self.output_dir, 'quant.genes.sf')
        }

        open(os.path.join(self.work_dir, '{}.here'.format(self.option('sample'))), 'w').close()

    def run(self):
        super(SalmonTool, self).run()
        self.run_salmon_quant()
        self.set_output()
        self.end()

    def run_salmon_quant(self):
        cmd = '{} quant'.format(self.program['salmon'])
        cmd += ' -l {}'.format(self.option('lib_type'))
        cmd += ' -i {}'.format(self.option('index').path)
        if self.option('unmated').is_set:
            cmd += ' -r {}'.format(self.option('unmated').path)
        else:
            cmd += ' -1 {}'.format(self.option('mates1').path)
            cmd += ' -2 {}'.format(self.option('mates2').path)
            if self.option('gc_bias'):
                cmd += ' --gcBias'
        cmd += ' -o {}'.format(self.output_dir)
        cmd += ' -p {}'.format(self.option('threads'))
        cmd += ' -g {}'.format(self.option('gene_map').path)
        if self.option('validate_mappings'):
            cmd += ' --validateMappings'
        if self.option('write_unmapped_names'):
            cmd += ' --writeUnmappedNames'
        runcmd(self, 'run_salmon_quant', cmd)

    def set_output(self):
        self.option('quant_sf').set_path(self.file['quant_sf'])
        self.option('quant_genes_sf').set_path(self.file['quant_genes_sf'])


class TestFunction(unittest.TestCase):
    '''
    This is test for the tool. Just run this script to do test.
    '''

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'salmon_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'tool',
            'name': 'whole_transcriptome.expression.salmon',
            'instant': False,
            'options': {
                'index': '/mnt/ilustre/users/sanger-dev/workspace/20190906/Single_index_4349_5691/Index/output/salmon.index',
                'mates1': '/mnt/ilustre/users/sanger-dev/workspace/20190905/WholeTranscriptome_workflow_3507_2079/FastpRna/output/fastq/S1.clean.1.fastq',
                'mates2': '/mnt/ilustre/users/sanger-dev/workspace/20190905/WholeTranscriptome_workflow_3507_2079/FastpRna/output/fastq/S1.clean.2.fastq',
                'threads': 16,
                'gene_map': '/mnt/ilustre/users/sanger-dev/workspace/20190906/Single_large_gush_5347_9539/LargeGush/MergeKnownNew/output/t2g.txt',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)

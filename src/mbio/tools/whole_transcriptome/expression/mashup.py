# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import os
import unittest

import pandas as pd
from biocluster.agent import Agent
from biocluster.tool import Tool

from mbio.packages.whole_transcriptome.utils import runcmd


class MashupAgent(Agent):
    '''
    last_modify: 2019.09.19
    '''

    def __init__(self, parent):
        super(MashupAgent, self).__init__(parent)
        options = [
            {'name': 'program', 'type': 'string', 'default': 'salmon'},
            {'name': 't_quant_list', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'g_quant_list', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'filter', 'type': 'bool', 'default': True},
            {'name': 'threshold', 'type': 'float', 'default': 0.0},
            {'name': 't2g', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 't_count', 'type': 'outfile', 'format': 'whole_transcriptome.common'},
            {'name': 'g_count', 'type': 'outfile', 'format': 'whole_transcriptome.common'},
        ]
        self.add_option(options)

    def check_options(self):
        for k, v in self._options.items():
            self.logger.debug('{} -> {}'.format(k, v.value))
        else:
            return True

    def set_resource(self):
        self._cpu = 1
        self._memory = '5G'

    def end(self):
        super(MashupAgent, self).end()


class MashupTool(Tool):
    def __init__(self, config):
        super(MashupTool, self).__init__(config)
        self.program = {
            'python': 'program/Python/bin/python'
        }
        self.script = {
            'mashup': os.path.join(self.config.PACKAGE_DIR, 'whole_transcriptome/expression/mashup.py')
        }
        self.file = {
            't_tpm': os.path.join(self.output_dir, 'T.tpm.txt'),
            'g_tpm': os.path.join(self.output_dir, 'G.tpm.txt'),
            't_fpkm': os.path.join(self.output_dir, 'T.fpkm.txt'),
            'g_fpkm': os.path.join(self.output_dir, 'G.fpkm.txt'),
            't_count': os.path.join(self.output_dir, 'T.count.txt'),
            'g_count': os.path.join(self.output_dir, 'G.count.txt')
        }

    def run(self):
        super(MashupTool, self).run()
        self.run_mashup_t()
        self.run_mashup_g()
        if self.option('filter'):
            self.exe_filter()
        self.set_output()
        self.end()

    def run_mashup_t(self):
        cmd = '{} {}'.format(self.program['python'], self.script['mashup'])
        cmd += ' -p {}'.format(self.option('program'))
        cmd += ' -i {}'.format(self.option('t_quant_list').path)
        cmd += ' -l {}'.format('T')
        cmd += ' -o {}'.format(self.output_dir)
        runcmd(self, 'run_mashup_t', cmd, block=False)

    def run_mashup_g(self):
        cmd = '{} {}'.format(self.program['python'], self.script['mashup'])
        cmd += ' -p {}'.format(self.option('program'))
        cmd += ' -i {}'.format(self.option('g_quant_list').path)
        cmd += ' -l {}'.format('G')
        cmd += ' -o {}'.format(self.output_dir)
        runcmd(self, 'run_mashup_g', cmd)

    def exe_filter(self):
        self.logger.info('start filtering expression matrix by threshold ({})'.format(self.option('threshold')))
        g_tpm_df = pd.read_table(self.file['g_tpm'], index_col=0)
        g_fpkm_df = pd.read_table(self.file['g_fpkm'], index_col=0)
        g_count_df = pd.read_table(self.file['g_count'], index_col=0)
        get_filtered_index = lambda g_df, threshold: g_df.loc[g_df.sum(axis=1) / g_df.shape[1] > threshold].index
        g_tpm_filtered_index, g_fpkm_filtered_index, g_count_filtered_index = \
            map(get_filtered_index, (g_tpm_df, g_fpkm_df, g_count_df), [self.option('threshold')] * 3)
        g_filtered_index = g_tpm_filtered_index & g_fpkm_filtered_index & g_count_filtered_index
        self.logger.info('succeed in obtaining {} gene id after screening'.format(len(g_filtered_index)))
        t_filtered_set = set()
        for line in open(self.option('t2g').path):
            t_id, g_id = line.strip().split('\t')
            if g_id in g_filtered_index:
                t_filtered_set.add(t_id)
        else:
            self.logger.info('succeed in obtaining {} transcript id after screening'.format(len(t_filtered_set)))
        for level in ['G', 'T']:
            for way in ['tpm', 'fpkm', 'count']:
                idx = {'G': g_filtered_index, 'T': t_filtered_set}[level]
                exp_df = pd.read_table(self.file['{}_{}'.format(level.lower(), way)], index_col=0)
                exp_df = exp_df.reindex(idx).fillna(float())
                exp_df.to_csv(self.file['{}_{}'.format(level.lower(), way)], sep='\t')
        else:
            self.logger.info('succeed in exporting filtered expression matrix')

    def set_output(self):
        self.option('t_count').set_path(self.file['t_count'])
        self.option('g_count').set_path(self.file['g_count'])


class TestFunction(unittest.TestCase):
    '''
    This is test for the tool. Just run this script to do test.
    '''

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'mashup_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'tool',
            'name': 'whole_transcriptome.expression.mashup',
            'instant': False,
            'options': {
                'program': 'salmon',
                't_quant_list': '/mnt/ilustre/users/sanger-dev/workspace/20190917/WholeTranscriptome_workflow_6312_7216/LargeGush/Expression/T.quant.list',
                'g_quant_list': '/mnt/ilustre/users/sanger-dev/workspace/20190917/WholeTranscriptome_workflow_6312_7216/LargeGush/Expression/G.quant.list',
                'filter': True,
                'threshold': 0.0,
                't2g': '/mnt/ilustre/users/sanger-dev/workspace/20190917/WholeTranscriptome_workflow_6312_7216/LargeGush/MergeKnownNew/output/t2g.txt',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)

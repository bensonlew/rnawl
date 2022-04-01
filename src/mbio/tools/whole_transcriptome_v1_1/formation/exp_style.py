# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import json
import os
import unittest

from biocluster.agent import Agent
from biocluster.tool import Tool

from mbio.packages.whole_transcriptome.utils import runcmd


class ExpStyleAgent(Agent):
    '''
    last_modify: 2019.01.09
    '''

    def __init__(self, parent):
        super(ExpStyleAgent, self).__init__(parent)
        options = [
            {'name': 't_tpm', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'g_tpm', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 't_fpkm', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'g_fpkm', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 't_style', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'g_style', 'type': 'infile', 'format': 'whole_transcriptome.common'},
        ]
        self.add_option(options)

    def check_options(self):
        for k, v in self._options.items():
            self.logger.debug('{} -> {}'.format(k, v.value))
        return True

    def set_resource(self):
        self._cpu = 1
        self._memory = '20G'

    def end(self):
        super(ExpStyleAgent, self).end()


class ExpStyleTool(Tool):
    def __init__(self, config):
        super(ExpStyleTool, self).__init__(config)
        self.program = {
            'python': 'program/Python/bin/python'
        }
        self.script = {
            'exp_style': os.path.join(self.config.PACKAGE_DIR, 'whole_transcriptome/formation/exp_style.py')
        }
        self.file = {
            'json': os.path.join('input.json')
        }

    def run(self):
        super(ExpStyleTool, self).run()
        self.pre_exp_style()
        self.run_exp_style()
        self.set_output()
        self.end()

    def pre_exp_style(self):
        json.dump({
            't_tpm': self.option('t_tpm').path,
            'g_tpm': self.option('g_tpm').path,
            't_fpkm': self.option('t_fpkm').path,
            'g_fpkm': self.option('g_fpkm').path,
            't_style': self.option('t_style').path,
            'g_style': self.option('g_style').path,
        }, open(self.file['json'], 'w'), indent=4)

    def run_exp_style(self):
        cmd = '{} {}'.format(self.program['python'], self.script['exp_style'])
        cmd += ' -i {}'.format(self.file['json'])
        cmd += ' -o {}'.format(self.output_dir)
        runcmd(self, 'run_exp_style', cmd)

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
            'id': 'exp_pca_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'tool',
            'name': 'whole_transcriptome.formation.exp_style',
            'instant': False,
            'options': {
                't_tpm': '/mnt/lustre/users/sanger/workspace/20200107/Longrna_sanger_233601/MessFlush/Expression/output/T.tpm.txt',
                'g_tpm': '/mnt/lustre/users/sanger/workspace/20200107/Longrna_sanger_233601/MessFlush/Expression/output/G.tpm.txt',
                't_fpkm': '/mnt/lustre/users/sanger/workspace/20200107/Longrna_sanger_233601/MessFlush/Expression/output/T.fpkm.txt',
                'g_fpkm': '/mnt/lustre/users/sanger/workspace/20200107/Longrna_sanger_233601/MessFlush/Expression/output/G.fpkm.txt',
                't_style': '/mnt/lustre/users/sanger/workspace/20200107/Longrna_sanger_233601/Assembly/output/trans_type.xls',
                'g_style': '/mnt/lustre/users/sanger/workspace/20200107/Longrna_sanger_233601/Assembly/output/gene_type.xls',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)

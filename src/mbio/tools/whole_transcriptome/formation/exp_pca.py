# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import os
import unittest

from biocluster.agent import Agent
from biocluster.tool import Tool

from mbio.packages.whole_transcriptome.utils import runcmd


class ExpPcaAgent(Agent):
    '''
    last_modify: 2019.10.22
    '''

    def __init__(self, parent):
        super(ExpPcaAgent, self).__init__(parent)
        options = [
            {'name': 'exp_matrix', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'group_table', 'type': 'infile', 'format': 'sample.group_table'},
            {'name': 'take_mean', 'type': 'bool', 'default': False},
            {'name': 'pca_table', 'type': 'outfile', 'format': 'whole_transcriptome.common'},
            {'name': 'evr_table', 'type': 'outfile', 'format': 'whole_transcriptome.common'},
            {'name': 'level', 'type': 'float', 'default': 0.95},
            {'name': 'ellipse_table', 'type': 'outfile', 'format': 'whole_transcriptome.common'},
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
        super(ExpPcaAgent, self).end()


class ExpPcaTool(Tool):
    def __init__(self, config):
        super(ExpPcaTool, self).__init__(config)
        self.program = {
            'python': 'program/Python/bin/python',
            'rscript': 'bioinfo/rconda/bin/Rscript',
        }
        self.script = {
            'exp_pca': os.path.join(self.config.PACKAGE_DIR, 'whole_transcriptome/formation/exp_pca.py'),
            'ellipse': os.path.join(self.config.PACKAGE_DIR, 'graph/scripts/Ellipse.R')
        }
        self.file = {
            'pca': os.path.join(self.output_dir, 'pca.txt'),
            'evr': os.path.join(self.output_dir, 'explained_variance_ratio.txt'),
            'ellipse': os.path.join(self.output_dir, 'ellipse.txt')
        }
        for group, samples in self.option('group_table').prop['group_dict'].items():
            if len(samples) < 3:
                self.flag_ellipse = False
                break
        else:
            self.flag_ellipse = True
        if len(self.option('group_table').prop['group_dict']) == 1:
            self.flag_ellipse = False

    def run(self):
        super(ExpPcaTool, self).run()
        self.run_exp_pca()
        if not self.option('take_mean') and self.flag_ellipse:
            self.run_ellipse()
        self.set_output()
        self.end()

    def run_exp_pca(self):
        cmd = '{} {}'.format(self.program['python'], self.script['exp_pca'])
        cmd += ' -i {}'.format(self.option('exp_matrix').path)
        cmd += ' -g {}'.format(self.option('group_table').path)
        if self.option('take_mean'):
            cmd += ' --mean'
        cmd += ' -o {}'.format(self.output_dir)
        runcmd(self, 'run_exp_pca', cmd)

    def run_ellipse(self):
        cmd = '{} {}'.format(self.program['rscript'], self.script['ellipse'])
        cmd += ' -f {}'.format(self.file['pca'])
        cmd += ' -l {}'.format(self.option('level'))
        cmd += ' -g {}'.format(self.option('group_table').path)
        cmd += ' -o {}'.format(self.file['ellipse'])
        runcmd(self, 'run_ellipse', cmd)

    def set_output(self):
        self.option('pca_table').set_path(self.file['pca'])
        self.option('evr_table').set_path(self.file['evr'])
        if not self.option('take_mean') and self.flag_ellipse:
            self.option('ellipse_table').set_path(self.file['ellipse'])


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
            'name': 'whole_transcriptome.formation.exp_pca',
            'instant': False,
            'options': {
                'exp_matrix': '/mnt/ilustre/users/sanger-dev/workspace/20191022/Single_exp_make_7164_9741/ExpMake/output/longrna/T.tpm.txt',
                'group_table': '/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/whole_transcriptome/demo_human/group.txt',
                'take_mean': False
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)

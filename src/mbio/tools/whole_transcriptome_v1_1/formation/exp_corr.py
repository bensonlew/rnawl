# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import os
import unittest

from biocluster.agent import Agent
from biocluster.tool import Tool

from mbio.packages.whole_transcriptome.utils import runcmd


class ExpCorrAgent(Agent):
    '''
    last_modify: 2019.09.26
    '''

    def __init__(self, parent):
        super(ExpCorrAgent, self).__init__(parent)
        CORR_METHOD = ('pearson', 'kendall', 'spearman')
        DIST_METHOD = ('euclidean', 'manhattan')
        CLUS_METHOD = ('complete', 'single', 'average')
        options = [
            {'name': 'exp_matrix', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'group_table', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'take_mean', 'type': 'bool', 'default': False},
            {'name': 'corr_method', 'type': 'string', 'default': CORR_METHOD[0]},
            {'name': 'dist_method', 'type': 'string', 'default': DIST_METHOD[0]},
            {'name': 'clus_method', 'type': 'string', 'default': CLUS_METHOD[0]},
            {'name': 'take_log', 'type': 'bool', 'default': False},
            {'name': 'corr_table', 'type': 'outfile', 'format': 'whole_transcriptome.common'},
            {'name': 'tree_file', 'type': 'outfile', 'format': 'whole_transcriptome.common'},
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
        super(ExpCorrAgent, self).end()


class ExpCorrTool(Tool):
    def __init__(self, config):
        super(ExpCorrTool, self).__init__(config)
        self.program = {
            'python': 'program/Python/bin/python',
            'rscript': 'bioinfo/rconda/bin/Rscript',
        }
        self.script = {
            'exp_corr': os.path.join(self.config.PACKAGE_DIR, 'whole_transcriptome/formation/exp_corr.py'),
            'hclustree': os.path.join(self.config.PACKAGE_DIR, 'whole_transcriptome/formation/hclustree.r'),
        }
        self.file = {
            'corr': os.path.join(self.output_dir, 'corr.txt'),
            'delim': os.path.join(self.output_dir, 'delim.txt'),
            'tree': os.path.join(self.output_dir, 'tree.txt')
        }

    def run(self):
        super(ExpCorrTool, self).run()
        self.run_exp_corr()
        self.run_hclustree()
        self.set_output()
        self.end()

    def run_exp_corr(self):
        cmd = '{} {}'.format(self.program['python'], self.script['exp_corr'])
        cmd += ' -i {}'.format(self.option('exp_matrix').path)
        cmd += ' -g {}'.format(self.option('group_table').path)
        cmd += ' -m {}'.format(self.option('corr_method'))
        if self.option('take_mean'):
            cmd += ' --mean'
        if self.option('corr_method') == 'pearson' and self.option('take_log'):
            cmd += ' --log'
        cmd += ' -o {}'.format(self.output_dir)
        runcmd(self, 'run_exp_corr', cmd)

    def run_hclustree(self):
        cmd = '{} {}'.format(self.program['rscript'], self.script['hclustree'])
        cmd += ' -i {}'.format(self.file['delim'])
        cmd += ' -d {}'.format(self.option('dist_method'))
        cmd += ' -c {}'.format(self.option('clus_method'))
        cmd += ' -o {}'.format(self.file['tree'])
        runcmd(self, 'run_hclustree', cmd)

    def set_output(self):
        self.option('corr_table').set_path(self.file['corr'])
        self.option('tree_file').set_path(self.file['tree'])


class TestFunction(unittest.TestCase):
    '''
    This is test for the tool. Just run this script to do test.
    '''

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'exp_corr_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'tool',
            'name': 'whole_transcriptome.formation.exp_corr',
            'instant': False,
            'options': {
                'exp_matrix': '/mnt/ilustre/users/sanger-dev/workspace/20190924/WholeTranscriptome_workflow_1623_1285/LargeGush/output/filter_by_express/classifyquant/T.tpm.txt',
                'group_table': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/whole_transcriptome/MJ20190118060/largeRNA/example_group.txt',
                'take_mean': False,
                'corr_method': 'pearson',
                'dist_method': 'euclidean',
                'clus_method': 'complete',
                'take_log': False
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)

# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import os
import shutil
import unittest

from biocluster.module import Module


class FormationModule(Module):
    def __init__(self, work_id):
        super(FormationModule, self).__init__(work_id)
        CORR_METHOD = ('pearson', 'kendall', 'spearman')
        DIST_METHOD = ('euclidean', 'manhattan')
        CLUS_METHOD = ('complete', 'single', 'average')
        options = [
            {'name': 'exp_matrix', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'group_table', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'threshold', 'type': 'float', 'default': 1.0},
            {'name': 'corr_method', 'type': 'string', 'default': CORR_METHOD[0]},
            {'name': 'dist_method', 'type': 'string', 'default': DIST_METHOD[0]},
            {'name': 'clus_method', 'type': 'string', 'default': CLUS_METHOD[0]},
            {'name': 'log', 'type': 'bool', 'default': False},
            {'name': 'take_mean', 'type': 'bool', 'default': False},
            {'name': 'confidence', 'type': 'float', 'default': 0.95},
        ]
        self.add_option(options)
        self.tools = list()

    def check_options(self):
        for k, v in self._options.items():
            self.logger.debug('{} -> {}'.format(k, v.value))
        else:
            return True

    def run(self):
        super(FormationModule, self).run()
        self.exp_graph = self.add_tool('whole_transcriptome.formation.exp_graph')
        self.exp_venn = self.add_tool('whole_transcriptome.formation.exp_venn')
        self.exp_corr = self.add_tool('whole_transcriptome.formation.exp_corr')
        self.exp_pca = self.add_tool('whole_transcriptome.formation.exp_pca')
        self.on_rely([self.exp_graph, self.exp_venn, self.exp_corr, self.exp_pca], self.set_output)
        self.run_exp_graph()
        self.run_exp_venn()
        self.run_exp_corr()
        self.run_exp_pca()

    def run_exp_graph(self):
        opts = {
            'exp_matrix': self.option('exp_matrix'),
            'group_table': self.option('group_table')
        }
        self.exp_graph.set_options(opts)
        self.exp_graph.run()

    def run_exp_venn(self):

        opts = {
            'exp_matrix': self.option('exp_matrix'),
            'group_table': self.option('group_table'),
            'threshold': self.option('threshold')
        }
        self.exp_venn.set_options(opts)
        self.exp_venn.run()

    def run_exp_corr(self):
        opts = {
            'exp_matrix': self.option('exp_matrix'),
            'group_table': self.option('group_table'),
            'corr_method': self.option('corr_method'),
            'dist_method': self.option('dist_method'),
            'clus_method': self.option('clus_method'),
            'log': self.option('log')
        }
        self.exp_corr.set_options(opts)
        self.exp_corr.run()

    def run_exp_pca(self):
        opts = {
            'exp_matrix': self.option('exp_matrix'),
            'group_table': self.option('group_table'),
            'take_mean': self.option('take_mean'),
            'level': self.option('confidence')
        }
        self.exp_pca.set_options(opts)
        self.exp_pca.run()

    def set_output(self):
        if os.path.isdir(self.output_dir):
            shutil.rmtree(self.output_dir)
        os.mkdir(self.output_dir)
        for name in ['exp_graph', 'exp_venn', 'exp_corr', 'exp_pca']:
            tool = getattr(self, name)
            shutil.copytree(tool.output_dir, os.path.join(self.output_dir, name))
        else:
            self.end()

    def end(self):
        super(FormationModule, self).end()


class TestFunction(unittest.TestCase):
    '''
    This is test for the module. Just run this script to do test.
    '''

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'formation_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'module',
            'name': 'whole_transcriptome.formation',
            'instant': False,
            'options': {
                'exp_matrix': '/mnt/ilustre/users/sanger-dev/workspace/20190924/Single_exp_make_9506_5592/ExpMake/output/mrna/G.tpm.txt',
                'group_table': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/whole_transcriptome/MJ20190118060/largeRNA/example_group.txt',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)

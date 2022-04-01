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
            {'name': 'calls', 'type': 'string', 'default': 'graph,venn,corr,pca'},
            {'name': 'graph_log', 'type': 'bool', 'default': True},
            {'name': 'threshold', 'type': 'float', 'default': 1.0},
            {'name': 'corr_method', 'type': 'string', 'default': CORR_METHOD[0]},
            {'name': 'dist_method', 'type': 'string', 'default': DIST_METHOD[0]},
            {'name': 'clus_method', 'type': 'string', 'default': CLUS_METHOD[0]},
            {'name': 'take_log', 'type': 'bool', 'default': False},
            {'name': 'take_mean', 'type': 'bool', 'default': False},
            {'name': 'confidence', 'type': 'float', 'default': 0.95},
        ]
        self.add_option(options)

    def check_options(self):
        for k, v in self._options.items():
            self.logger.debug('{} -> {}'.format(k, v.value))
        else:
            return True

    def run(self):
        super(FormationModule, self).run()
        self.tools = dict()
        self.funcs = list()
        if 'graph' in self.option('calls'):
            self.tools['exp_graph'] = self.add_tool('whole_transcriptome.formation.exp_graph')
            self.funcs.append(self.run_exp_graph)
        if 'venn' in self.option('calls'):
            self.tools['exp_venn'] = self.add_tool('whole_transcriptome.formation.exp_venn')
            self.funcs.append(self.run_exp_venn)
        if 'corr' in self.option('calls'):
            self.tools['exp_corr'] = self.add_tool('whole_transcriptome.formation.exp_corr')
            self.funcs.append(self.run_exp_corr)
        if 'pca' in self.option('calls'):
            self.tools['exp_pca'] = self.add_tool('whole_transcriptome.formation.exp_pca')
            self.funcs.append(self.run_exp_pca)
        self.on_rely(self.tools.values(), self.set_output)
        for func in self.funcs:
            func()

    def run_exp_graph(self):
        opts = {
            'exp_matrix': self.option('exp_matrix'),
            'group_table': self.option('group_table')
        }
        self.tools['exp_graph'].set_options(opts)
        self.tools['exp_graph'].run()

    def run_exp_venn(self):
        opts = {
            'exp_matrix': self.option('exp_matrix'),
            'group_table': self.option('group_table'),
            'threshold': self.option('threshold')
        }
        self.tools['exp_venn'].set_options(opts)
        self.tools['exp_venn'].run()

    def run_exp_corr(self):
        opts = {
            'exp_matrix': self.option('exp_matrix'),
            'group_table': self.option('group_table'),
            'corr_method': self.option('corr_method'),
            'dist_method': self.option('dist_method'),
            'clus_method': self.option('clus_method'),
            'take_log': self.option('take_log')
        }
        self.tools['exp_corr'].set_options(opts)
        self.tools['exp_corr'].run()

    def run_exp_pca(self):
        opts = {
            'exp_matrix': self.option('exp_matrix'),
            'group_table': self.option('group_table').path,
            'take_mean': self.option('take_mean'),
            'level': self.option('confidence')
        }
        self.tools['exp_pca'].set_options(opts)
        self.tools['exp_pca'].run()

    def set_output(self):
        if os.path.isdir(self.output_dir):
            shutil.rmtree(self.output_dir)
        os.mkdir(self.output_dir)
        for name, tool in self.tools.items():
            shutil.copytree(tool.output_dir, os.path.join(self.output_dir, name))
        self.end()

    def end(self):
        super(FormationModule, self).end()


class TestFunction(unittest.TestCase):
    '''
    This is test for the module. Just run this script to do test.
    '''

    def test_mrna(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'formation_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'module',
            'name': 'whole_transcriptome.formation',
            'instant': False,
            'options': {
                'exp_matrix': '/mnt/ilustre/users/sanger-dev/workspace/20191106/WholeTranscriptome_tsg_36088/Transfer/output/exp_make/mrna/T.tpm.txt',
                'group_table': '/mnt/ilustre/users/sanger-dev/workspace/20191106/WholeTranscriptome_tsg_36088/Transfer/output/group.txt',
                'calls': 'graph,venn'
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

    def test_lncrna(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'formation_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'module',
            'name': 'whole_transcriptome.formation',
            'instant': False,
            'options': {
                'exp_matrix': '/mnt/ilustre/users/sanger-dev/workspace/20191106/WholeTranscriptome_tsg_36088/Transfer/output/exp_make/lncrna/T.tpm.txt',
                'group_table': '/mnt/ilustre/users/sanger-dev/workspace/20191106/WholeTranscriptome_tsg_36088/Transfer/output/group.txt',
                'calls': 'graph,venn'
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

    def test_circrna(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'formation_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'module',
            'name': 'whole_transcriptome.formation',
            'instant': False,
            'options': {
                'exp_matrix': '/mnt/ilustre/users/sanger-dev/workspace/20191106/WholeTranscriptome_tsg_36088/Transfer/output/exp_make/circrna/T.rpm.txt',
                'group_table': '/mnt/ilustre/users/sanger-dev/workspace/20191106/WholeTranscriptome_tsg_36088/Transfer/output/group.txt',
                'calls': 'graph,venn'
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

    def test_mirna(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'formation_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'module',
            'name': 'whole_transcriptome.formation',
            'instant': False,
            'options': {
                'exp_matrix': '/mnt/ilustre/users/sanger-dev/workspace/20191106/WholeTranscriptome_tsg_36088/Transfer/output/exp_make/mirna/T.tpm.txt',
                'group_table': '/mnt/ilustre/users/sanger-dev/workspace/20191106/WholeTranscriptome_tsg_36088/Transfer1/output/group.txt',
                'calls': 'graph,venn'
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test_mirna')])
    unittest.TextTestRunner(verbosity=2).run(suite)

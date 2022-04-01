# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import os
import unittest

from biocluster.module import Module


class DiffExpModule(Module):
    def __init__(self, work_id):
        super(DiffExpModule, self).__init__(work_id)
        PROGRAM = ('DESeq2', 'edgeR', 'DEGseq')
        FILTER = ('none', 'mean', 'max', 'min', 'sum', 'median')
        METHOD = ('bonferroni', 'holm', 'bh', 'by')
        STAT_TYPE = ('pvalue', 'padjust')
        options = [
            {'name': 'program', 'type': 'string', 'default': PROGRAM[0]},
            {'name': 'count_matrix', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'group_table', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'control_table', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'exp_matrix', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'kind_table', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'filter', 'type': 'string', 'default': FILTER[0]},
            {'name': 'threshold', 'type': 'float', 'default': 0.0},
            {'name': 'method', 'type': 'string', 'default': METHOD[2]},
            {'name': 'stat_type', 'type': 'string', 'default': STAT_TYPE[1]},
            {'name': 'stat_cutoff', 'type': 'float', 'default': 0.05},
            {'name': 'fc', 'type': 'float', 'default': 2.0},
            {'name': 'email', 'type': 'bool', 'default': False},
        ]
        self.add_option(options)
        self.tools = list()

    def check_options(self):
        for k, v in self._options.items():
            self.logger.debug('{} -> {}'.format(k, v.value))
        return True

    def run(self):
        super(DiffExpModule, self).run()
        if self.option('program') == 'DESeq2':
            self.run_deseq2()
        elif self.option('program') == 'edgeR':
            self.run_edger()
        elif self.option('program') == 'DEGseq':
            self.run_degseq()

    def run_deseq2(self):
        self.deseq2 = self.add_tool('whole_transcriptome.expression.deseq2')
        opts = {
            'count_matrix': self.option('count_matrix'),
            'group_table': self.option('group_table'),
            'control_table': self.option('control_table')
        }
        self.deseq2.set_options(opts)
        self.deseq2.on('end', self.run_uniform)
        self.deseq2.run()

    def run_edger(self):
        self.edger = self.add_tool('whole_transcriptome.expression.edger')
        opts = {
            'count_matrix': self.option('count_matrix'),
            'group_table': self.option('group_table'),
            'control_table': self.option('control_table')
        }
        self.edger.set_options(opts)
        self.edger.on('end', self.run_uniform)
        self.edger.run()

    def run_degseq(self):
        self.degseq = self.add_tool('whole_transcriptome.expression.degseq')
        opts = {
            'count_matrix': self.option('count_matrix'),
            'group_table': self.option('group_table'),
            'control_table': self.option('control_table')
        }
        self.degseq.set_options(opts)
        self.degseq.on('end', self.run_uniform)
        self.degseq.run()

    def run_uniform(self):
        self.uniform = self.add_tool('whole_transcriptome.expression.uniform')
        if self.option('program') == 'DESeq2':
            input_dir = self.deseq2.output_dir
        elif self.option('program') == 'edgeR':
            input_dir = self.edger.output_dir
        elif self.option('program') == 'DEGseq':
            input_dir = self.degseq.output_dir
        opts = {
            'program': self.option('program'),
            'input_dir': input_dir,
            'group_table': self.option('group_table'),
            'control_table': self.option('control_table'),
            'exp_matrix': self.option('exp_matrix'),
            'kind_table': self.option('kind_table'),
            'filter': self.option('filter'),
            'threshold': self.option('threshold'),
            'method': self.option('method'),
            'stat_type': self.option('stat_type'),
            'stat_cutoff': self.option('stat_cutoff'),
            'fc': self.option('fc')
        }
        self.uniform.set_options(opts)
        self.uniform.on('end', self.set_output)
        self.uniform.run()

    def set_output(self):
        for file_name in os.listdir(self.uniform.output_dir):
            source = os.path.join(self.uniform.output_dir, file_name)
            link_name = os.path.join(self.output_dir, file_name)
            if os.path.isfile(link_name):
                os.remove(link_name)
            os.link(source, link_name)
        self.option('email', self.uniform.option('email'))
        self.end()

    def end(self):
        super(DiffExpModule, self).end()


class TestFunction(unittest.TestCase):
    '''
    This is test for the module. Just run this script to do test.
    '''

    def test_longrna(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'diff_exp_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'module',
            'name': 'whole_transcriptome.diff_exp',
            'instant': False,
            'options': {
                'count_matrix': '/mnt/ilustre/users/sanger-dev/workspace/20191031/Longrna_tsg_36023/ExpMake/output/count/T.reads.txt',
                'group_table': '/mnt/ilustre/users/sanger-dev/workspace/20191031/Longrna_tsg_36023/remote_input/group_table/group.txt',
                'control_table': '/mnt/ilustre/users/sanger-dev/workspace/20191031/Longrna_tsg_36023/remote_input/control_file/control.txt',
                'exp_matrix': '/mnt/ilustre/users/sanger-dev/workspace/20191031/Longrna_tsg_36023/ExpMake/output/longrna/T.tpm.txt',
                'kind_table': '/mnt/ilustre/users/sanger-dev/workspace/20191031/Longrna_tsg_36023/DiffExp/kind.txt',
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
            'id': 'diff_exp_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'module',
            'name': 'whole_transcriptome.diff_exp',
            'instant': False,
            'options': {
                'count_matrix': '/mnt/ilustre/users/sanger-dev/workspace/20191031/Longrna_tsg_36023/ExpMake/output/count/C.reads.txt',
                'group_table': '/mnt/ilustre/users/sanger-dev/workspace/20191031/Longrna_tsg_36023/remote_input/group_table/group.txt',
                'control_table': '/mnt/ilustre/users/sanger-dev/workspace/20191031/Longrna_tsg_36023/remote_input/control_file/control.txt',
                'exp_matrix': '/mnt/ilustre/users/sanger-dev/workspace/20191031/Longrna_tsg_36023/ExpMake/output/circrna/T.rpm.txt',
                'kind_table': '/mnt/ilustre/users/sanger-dev/workspace/20191031/Longrna_tsg_36023/DiffExp1/kind.txt',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

    def test_edger(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'diff_exp_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'module',
            'name': 'whole_transcriptome.diff_exp',
            'instant': False,
            'options': {
                'program': 'edgeR',
                'count_matrix': '/mnt/ilustre/users/sanger-dev/workspace/20191031/Longrna_tsg_36023/ExpMake/output/count/T.reads.txt',
                'group_table': '/mnt/ilustre/users/sanger-dev/workspace/20191031/Longrna_tsg_36023/remote_input/group_table/group.txt',
                'control_table': '/mnt/ilustre/users/sanger-dev/workspace/20191031/Longrna_tsg_36023/remote_input/control_file/control.txt',
                'exp_matrix': '/mnt/ilustre/users/sanger-dev/workspace/20191031/Longrna_tsg_36023/ExpMake/output/longrna/T.tpm.txt',
                'kind_table': '/mnt/ilustre/users/sanger-dev/workspace/20191031/Longrna_tsg_36023/DiffExp/kind.txt',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

    def test_degseq(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'diff_exp_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'module',
            'name': 'whole_transcriptome.diff_exp',
            'instant': False,
            'options': {
                'program': 'DEGseq',
                'count_matrix': '/mnt/ilustre/users/sanger-dev/workspace/20191031/Longrna_tsg_36023/ExpMake/output/count/T.reads.txt',
                'group_table': '/mnt/ilustre/users/sanger-dev/workspace/20191031/Longrna_tsg_36023/remote_input/group_table/group.txt',
                'control_table': '/mnt/ilustre/users/sanger-dev/workspace/20191031/Longrna_tsg_36023/remote_input/control_file/control.txt',
                'exp_matrix': '/mnt/ilustre/users/sanger-dev/workspace/20191031/Longrna_tsg_36023/ExpMake/output/longrna/T.tpm.txt',
                'kind_table': '/mnt/ilustre/users/sanger-dev/workspace/20191031/Longrna_tsg_36023/DiffExp/kind.txt',
                'stat_cutoff': '0.001',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test_degseq')])
    unittest.TextTestRunner(verbosity=2).run(suite)

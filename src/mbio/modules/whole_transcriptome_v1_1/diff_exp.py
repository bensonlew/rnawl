# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import os
import unittest

from biocluster.module import Module


class DiffExpModule(Module):
    def __init__(self, work_id):
        super(DiffExpModule, self).__init__(work_id)
        PROGRAM = ('DESeq2', 'edgeR', 'DEGseq', 'limma', 'Noiseq')
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
            {'name': 'prob', 'type': 'float', 'default': 0.8}, #新增 NOIseq软件专用参数
            {'name': 'is_batch', 'type': 'bool', 'default': False},
            {'name': 'has_batch', 'type': 'bool', 'default': True},
            {'name': 'batch_matrix', 'type': 'infile', 'format': "whole_transcriptome.common"},
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
        elif self.option('program') == 'limma':
            self.run_limma()
        elif self.option('program').lower() == 'noiseq':
            self.run_noiseq()
        elif self.option('program').lower() == 'svaseqlimma':
            self.run_svalimma()
    def run_deseq2(self):
        self.deseq2 = self.add_tool('whole_transcriptome_v1_1.expression.deseq2')
        opts = {
            'count_matrix': self.option('count_matrix'),
            'group_table': self.option('group_table'),
            'control_table': self.option('control_table'),
            'is_batch': self.option('is_batch')
        }
        if self.option('is_batch') == False:
            pass
        if self.option('is_batch') == True and self.option('has_batch') == True:
            opts.update(
                has_batch=self.option('has_batch'),
                batch_matrix=self.option('batch_matrix')
            )
        self.deseq2.set_options(opts)
        self.deseq2.on('end', self.run_uniform)
        self.deseq2.run()

    def run_edger(self):
        self.edger = self.add_tool('whole_transcriptome_v1_1.expression.edger')
        opts = {
            'count_matrix': self.option('count_matrix'),
            'group_table': self.option('group_table'),
            'control_table': self.option('control_table'),
            'is_batch': self.option('is_batch')
        }
        if self.option('is_batch') == False:
            pass
        if self.option('is_batch') == True and self.option('has_batch') == True:
            opts.update(
                has_batch=self.option('has_batch'),
                batch_matrix=self.option('batch_matrix')
            )
        self.edger.set_options(opts)
        self.edger.on('end', self.run_uniform)
        self.edger.run()

    def run_degseq(self):
        self.degseq = self.add_tool('whole_transcriptome_v1_1.expression.degseq')
        opts = {
            'count_matrix': self.option('count_matrix'),
            'group_table': self.option('group_table'),
            'control_table': self.option('control_table')
        }
        self.degseq.set_options(opts)
        self.degseq.on('end', self.run_uniform)
        self.degseq.run()

    def run_limma(self):
        self.limma = self.add_tool('whole_transcriptome_v1_1.expression.limma')
        opts = {
            'count_matrix': self.option('count_matrix'),
            'group_table': self.option('group_table'),
            'control_table': self.option('control_table'),
            'is_batch': self.option('is_batch')
        }
        if self.option('is_batch') == False:
            pass
        if self.option('is_batch') == True and self.option('has_batch') == True:
            opts.update(
                has_batch=self.option('has_batch'),
                batch_matrix=self.option('batch_matrix')
            )
        self.limma.set_options(opts)
        self.limma.on('end', self.run_uniform)
        self.limma.run()

    def run_noiseq(self):
        self.noiseq = self.add_tool('whole_transcriptome_v1_1.expression.noiseq')
        opts = {
            'count_matrix': self.option('count_matrix'),
            'group_table': self.option('group_table'),
            'control_table': self.option('control_table')
        }
        self.noiseq.set_options(opts)
        self.noiseq.on('end', self.run_uniform)
        self.noiseq.run()

    def run_svalimma(self):
        self.svalimma = self.add_tool('whole_transcriptome_v1_1.expression.svalimma')
        opts = {
            'count_matrix': self.option('count_matrix'),
            'group_table': self.option('group_table'),
            'control_table': self.option('control_table')
        }
        self.svalimma.set_options(opts)
        self.svalimma.on('end', self.run_uniform)
        self.svalimma.run()

    def run_uniform(self):
        self.uniform = self.add_tool('whole_transcriptome_v1_1.expression.uniform')
        if self.option('program') == 'DESeq2':
            input_dir = self.deseq2.output_dir
        elif self.option('program') == 'edgeR':
            input_dir = self.edger.output_dir
        elif self.option('program') == 'DEGseq':
            input_dir = self.degseq.output_dir
        elif self.option('program') == 'limma':
            input_dir = self.limma.output_dir
        elif self.option('program').lower() == 'noiseq':
            input_dir = self.noiseq.output_dir
        elif self.option('program').lower() == 'svaseqlimma':
            input_dir = self.svalimma.output_dir
        opts = {
            'program': self.option('program'),
            'input_dir': input_dir,
            'group_table': self.option('group_table'),
            'control_table': self.option('control_table'),
            'exp_matrix': self.option('exp_matrix'),
            'kind_table': self.option('kind_table'),
            'filter': self.option('filter'),
            'threshold': self.option('threshold'),
            'fc': self.option('fc')
        }
        if self.option('program').lower() in ['deseq2', 'edger2', 'degseq', 'limma', 'svalimma']:
            opts.update({
                'method': self.option('method'),
                'stat_type': self.option('stat_type'),
                'stat_cutoff': self.option('stat_cutoff'),
            })
        elif self.option('program').lower() in ['noiseq']:
            opts.update({
                'prob': self.option('prob')
            })
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
            'name': 'whole_transcriptome2.diff_exp',
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

    def test_limma(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'diff_exp_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'module',
            'name': 'whole_transcriptome2.diff_exp',
            'instant': False,
            'options': {
                'program': 'limma',
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

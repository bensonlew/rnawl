# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import os
import unittest

import pandas as pd
from biocluster.agent import Agent
from biocluster.tool import Tool

from mbio.packages.whole_transcriptome.utils import runcmd


class UniformAgent(Agent):
    """
    last_modify: 2019.10.12
    """

    def __init__(self, parent):
        super(UniformAgent, self).__init__(parent)
        PROGRAM = ('DESeq2', 'edgeR', 'DEGseq', 'limma', 'NOIseq', 'svaseqlimma')
        FILTER = ('none', 'mean', 'max', 'min', 'sum', 'median')
        METHOD = ('bonferroni', 'holm', 'bh', 'by')
        STAT_TYPE = ('pvalue', 'padjust')
        options = [
            {'name': 'program', 'type': 'string', 'default': PROGRAM[0]},
            {'name': 'input_dir', 'type': 'infile', 'format': 'whole_transcriptome.common_dir'},
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
            {'name': 'prob', 'type': 'float', 'default': 0.8},
            {'name': 'email', 'type': 'bool', 'default': False},
        ]
        self.add_option(options)

    def check_options(self):
        for k, v in self._options.items():
            self.logger.debug('{} -> {}'.format(k, v.value))
        else:
            return True

    def set_resource(self):
        self._cpu = 1
        self._memory = '10G'

    def end(self):
        super(UniformAgent, self).end()


class UniformTool(Tool):
    def __init__(self, config):
        super(UniformTool, self).__init__(config)
        self.program = {
            'python': 'miniconda2/bin/python',
        }
        self.script = {
            'uniform': os.path.join(self.config.PACKAGE_DIR, 'whole_transcriptome_v1_1/expression/uniform.py')
        }

    def run(self):
        super(UniformTool, self).run()
        self.run_uniform()
        self.set_output()
        self.end()

    def run_uniform(self):
        cmd = '{} {}'.format(self.program['python'], self.script['uniform'])
        cmd += ' -u {}'.format(self.option('program'))
        cmd += ' -i {}'.format(self.option('input_dir').path)
        cmd += ' -e {}'.format(self.option('exp_matrix').path)
        cmd += ' -g {}'.format(self.option('group_table').path)
        cmd += ' -c {}'.format(self.option('control_table').path)
        cmd += ' -k {}'.format(self.option('kind_table').path)
        cmd += ' -f {}'.format(self.option('filter'))
        cmd += ' -t {}'.format(self.option('threshold'))
        cmd += ' -m {}'.format(self.option('method'))
        cmd += ' -s {}'.format(self.option('stat_type'))
        cmd += ' -p {}'.format(self.option('stat_cutoff'))
        cmd += ' -d {}'.format(self.option('fc'))
        cmd += ' -o {}'.format(self.output_dir)
        cmd += ' -prob {}'.format(self.option('prob'))
        runcmd(self, 'run_uniform', cmd)

    def set_output(self):
        df = pd.read_table(os.path.join(self.output_dir, 'summary.txt'))
        _list = list()
        for vs_pair in df.columns[1:-1]:
            n_yes = 0
            for v in df[vs_pair]:
                eles = str(v).split('|')
                if len(eles) == 2 and eles[0] == 'yes':
                    n_yes += 1
            _list.append(n_yes < 10)
        self.option('email', all(_list))


class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'uniform_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'tool',
            'name': 'whole_transcriptome.expression.uniform',
            'instant': False,
            'options': {
                'input_dir': '/mnt/ilustre/users/sanger-dev/workspace/20191031/Longrna_tsg_36023/DiffExp/Deseq2/output',
                'group_table': '/mnt/ilustre/users/sanger-dev/workspace/20191031/Longrna_tsg_36023/remote_input/group_table/group.txt',
                'control_table': '/mnt/ilustre/users/sanger-dev/workspace/20191031/Longrna_tsg_36023/remote_input/control_file/control.txt',
                'exp_matrix': '/mnt/ilustre/users/sanger-dev/workspace/20191031/Longrna_tsg_36023/ExpMake/output/longrna/T.tpm.txt',
                'kind_table': '/mnt/ilustre/users/sanger-dev/workspace/20191031/Longrna_tsg_36023/DiffExp/kind.txt',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)

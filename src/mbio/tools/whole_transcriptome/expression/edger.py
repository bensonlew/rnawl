# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import json
import os
import shutil
import unittest

import numpy as np
import pandas as pd
from biocluster.agent import Agent
from biocluster.tool import Tool

from mbio.packages.whole_transcriptome.utils import runcmd


class EdgerAgent(Agent):
    '''
    last_modify: 2019.11.15
    '''

    def __init__(self, parent):
        super(EdgerAgent, self).__init__(parent)
        FIT_METHOD = ('exactTest', 'glmFit', 'glmQLFit')
        options = [
            {'name': 'count_matrix', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'group_table', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'control_table', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'fit_method', 'type': 'string', 'default': FIT_METHOD[2]},
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
        super(EdgerAgent, self).end()


class EdgerTool(Tool):
    def __init__(self, config):
        super(EdgerTool, self).__init__(config)
        self.program = {
            'rscript': 'bioinfo/rconda/bin/Rscript',
        }
        self.script = {
            'edger': os.path.join(self.config.PACKAGE_DIR, 'whole_transcriptome/expression/edger.r')
        }
        self.file = {
            'json': os.path.join(self.work_dir, 'input.json'),
            'count': os.path.join(self.work_dir, 'count.txt'),
            'group': os.path.join(self.work_dir, 'group.txt'),
            'control': os.path.join(self.work_dir, 'control.txt')
        }

    def run(self):
        super(EdgerTool, self).run()
        self.pre_edger()
        self.run_edger()
        self.set_output()
        self.end()

    def pre_edger(self):
        shutil.copy(self.option('group_table').path, self.file['group'])
        shutil.copy(self.option('control_table').path, self.file['control'])
        samples = list()
        groups = list()
        for line in open(self.file['group']):
            if line.strip() and line[0] != '#':
                sample, group = line.strip().split('\t')
                if sample not in samples:
                    samples.append(sample)
                    groups.append(group)
        contrasts = list()
        for line in open(self.file['control']):
            if line.strip() and line[0] != '#':
                ctrl, case = line.strip().split('\t')
                contrasts.append({'ctrl': ctrl, 'case': case})
        df = pd.read_table(self.option('count_matrix').path, index_col=0)
        df = df.reindex(samples, axis=1)
        df.to_csv(self.file['count'], sep='\t')
        if len(set(samples)) == len(set(groups)):
            fit_method = 'exactTest'
        else:
            fit_method = self.option('fit_method')
        json.dump({
            'count': self.file['count'],
            'groups': groups,
            'contrasts': contrasts,
            'fit_method': fit_method,
            'output': self.output_dir
        }, open(self.file['json'], 'w'), indent=4)

    def run_edger(self):
        cmd = '{} {} -i {}'.format(self.program['rscript'], self.script['edger'], self.file['json'])
        runcmd(self, 'run_edger', cmd)

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
            'id': 'edger_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'tool',
            'name': 'whole_transcriptome.expression.edger',
            'instant': False,
            'options': {
                'count_matrix': '/mnt/ilustre/users/sanger-dev/workspace/20191114/DiffExp_tsg_36088_6566_5962/DiffPrepare/output/G.mRNA.all.txt',
                'group_table': '/mnt/ilustre/users/sanger-dev/workspace/20191114/DiffExp_tsg_36088_6566_5962/ExpBuild/output/group.txt',
                'control_table': '/mnt/ilustre/users/sanger-dev/workspace/20191114/DiffExp_tsg_36088_6566_5962/ExpBuild/output/control.txt'
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)

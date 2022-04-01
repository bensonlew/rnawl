# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import json
import os
import shutil
import unittest

import pandas as pd
from biocluster.agent import Agent
from biocluster.tool import Tool

from mbio.packages.whole_transcriptome.utils import runcmd


class DegseqAgent(Agent):
    '''
    last_modify: 2019.11.15
    '''

    def __init__(self, parent):
        super(DegseqAgent, self).__init__(parent)
        options = [
            {'name': 'count_matrix', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'group_table', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'control_table', 'type': 'infile', 'format': 'whole_transcriptome.common'},
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
        super(DegseqAgent, self).end()


class DegseqTool(Tool):
    def __init__(self, config):
        super(DegseqTool, self).__init__(config)
        self.program = {
            'rscript': 'bioinfo/rconda/bin/Rscript',
        }
        self.script = {
            'degseq': os.path.join(self.config.PACKAGE_DIR, 'whole_transcriptome_v1_1/expression/degseq.r')
        }
        self.file = {
            'json': os.path.join(self.work_dir, 'input.json'),
            'count': os.path.join(self.work_dir, 'count.txt'),
            'group': os.path.join(self.work_dir, 'group.txt'),
            'control': os.path.join(self.work_dir, 'control.txt')
        }

    def run(self):
        super(DegseqTool, self).run()
        self.pre_degseq()
        self.run_degseq()
        self.set_output()
        self.end()

    def pre_degseq(self):
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
        sample_column_dict = dict()
        for i, (sample, group) in enumerate(zip(samples, groups), start=2):
            if group in sample_column_dict:
                sample_column_dict[group].append(i)
            else:
                sample_column_dict[group] = [i]
        contrasts = list()
        columns = list()
        for line in open(self.file['control']):
            if line.strip() and line[0] != '#':
                ctrl, case = line.strip().split('\t')
                contrasts.append({'ctrl': ctrl, 'case': case})
                columns.append({'ctrl': sample_column_dict[ctrl], 'case': sample_column_dict[case]})
        df = pd.read_table(self.option('count_matrix').path, index_col=0)
        df = df.reindex(samples, axis=1)
        df.to_csv(self.file['count'], sep='\t')
        json.dump({
            'count': self.file['count'],
            'columns': columns,
            'contrasts': contrasts,
            'output': self.output_dir
        }, open(self.file['json'], 'w'), indent=4)

    def run_degseq(self):
        cmd = '{} {} -i {}'.format(self.program['rscript'], self.script['degseq'], self.file['json'])
        runcmd(self, 'run_degseq', cmd)

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
            'id': 'degseq_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'tool',
            'name': 'whole_transcriptome.expression.degseq',
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

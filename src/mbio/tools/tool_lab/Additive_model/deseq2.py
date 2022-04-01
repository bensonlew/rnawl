# -*- coding: utf-8 -*-
# __author__ = 'zoujiaxun'

import json
import os
import shutil
import unittest

import pandas as pd
from biocluster.agent import Agent
from biocluster.tool import Tool

from mbio.packages.whole_transcriptome.utils import runcmd


class Deseq2Agent(Agent):
    '''
    last_modify: 2019.12.12
    '''

    def __init__(self, parent):
        super(Deseq2Agent, self).__init__(parent)
        options = [
            {'name': 'batch', 'type': 'bool', 'default': False},
            {'name': 'paired', 'type': 'bool', 'default': False},
            {'name': 'count_matrix', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'group_table', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'control_table', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'together', 'type': 'bool', 'default': True},
            {'name': 'batch_matrix', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'paired_matrix', 'type': 'infile', 'format': 'whole_transcriptome.common'},
        ]
        self.add_option(options)

    def check_options(self):
        for k, v in self._options.items():
            self.logger.debug('{} -> {}'.format(k, v.value))
        else:
            return True

    def set_resource(self):
        self._cpu = 1
        self._memory = '20G'

    def end(self):
        super(Deseq2Agent, self).end()


class Deseq2Tool(Tool):
    def __init__(self, config):
        super(Deseq2Tool, self).__init__(config)
        self.program = {
            'python': 'program/Python/bin/python',
            'rscript': 'bioinfo/rconda/bin/Rscript',
        }
        self.script = {
            'deseq2': os.path.join(self.config.PACKAGE_DIR, 'tool_lab/Additive_model/deseq2.r')
        }
        self.file = {
            'json': os.path.join(self.work_dir, 'input.json'),
            'count': os.path.join(self.work_dir, 'count.txt'),
            'design': os.path.join(self.work_dir, 'design.txt'),
            'group': os.path.join(self.work_dir, 'group.txt'),
            'contrast': os.path.join(self.work_dir, 'contrast.txt'),
            'batch': os.path.join(self.work_dir, 'batch.txt'),
            'paired': os.path.join(self.work_dir, 'paired.txt')
        }

    def run(self):
        super(Deseq2Tool, self).run()
        self.pre_deseq2()
        self.run_deseq2()
        self.set_output()
        self.end()

    def pre_deseq2(self):
        shutil.copy(self.option('group_table').path, self.file['group'])
        shutil.copy(self.option('control_table').path, self.file['contrast'])
        batchs = list()
        samples = list()
        groups = list()
        paireds = list()
        if self.option('batch'):
            shutil.copy(self.option('batch_matrix').path, self.file['batch'])
            for line in open(self.file['batch']):
                if line.strip() and line[0] != '#':
                    sample, batch = line.strip().split('\t')
                    batchs.append(batch)
        if self.option('paired'):
            shutil.copy(self.option('paired_matrix').path, self.file['paired'])
            for line in open(self.file['paired']):
                if line.strip() and line[0] != '#':
                    sample, paired = line.strip().split('\t')
                    paireds.append(paired)

        for line in open(self.file['group']):
            if line.strip() and line[0] != '#':
                sample, group = line.strip().split('\t')
                if sample not in samples:
                    samples.append(sample)
                    groups.append(group)

        contrasts = list()
        for line in open(self.file['contrast']):
            if line.strip() and line[0] != '#':
                ctrl, case = line.strip().split('\t')
                contrasts.append({'ctrl': ctrl, 'case': case})
        df = pd.read_table(self.option('count_matrix').path, index_col=0)
        df = df.reindex(samples, axis=1)
        df = df.astype(int)
        df.to_csv(self.file['count'], sep='\t')
        json.dump({
            'Batch': self.option('batch'),
            'Paired': self.option('paired'),
            'count': self.file['count'],
            'samples': samples,
            'groups': groups,
            'batchs': batchs,
            'paireds': paireds,
            'contrasts': contrasts,
            'together': self.option('together'),
            'output': self.output_dir
        }, open(self.file['json'], 'w'), indent=4)

    def run_deseq2(self):
        cmd = '{} {} -i {}'.format(self.program['rscript'], self.script['deseq2'], self.file['json'])
        runcmd(self, 'run_deseq2', cmd)

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
            'id': 'deseq2_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'tool',
            'name': 'tool_lab.Additive_model.deseq2',
            'instant': False,
            'options': {
                'count_matrix': '/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/batch/gene.count.matrix.txt',
                'group_table': '/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/batch/group.txt',
                'control_table': '/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/batch/contrast.txt',
                'batch_matrix':'/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/batch/batch.txt',
                'together': 'True',
                'batch': True
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)

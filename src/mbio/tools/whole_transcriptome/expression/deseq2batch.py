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


class Deseq2batchAgent(Agent):
    '''
    last_modify: 2019.12.12
    '''

    def __init__(self, parent):
        super(Deseq2batchAgent, self).__init__(parent)
        options = [
            {'name': 'count_matrix', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'group_table', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'control_table', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'together', 'type': 'bool', 'default': True},
            {'name': 'pheno_matrix', 'type': 'infile', 'format': 'whole_transcriptome.common'}
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
        super(Deseq2batchAgent, self).end()


class Deseq2batchTool(Tool):
    def __init__(self, config):
        super(Deseq2batchTool, self).__init__(config)
        self.program = {
            'python': 'miniconda2/bin/python',
            'rscript': 'bioinfo/rconda/bin/Rscript',
        }
        self.script = {
            'deseq2batch': os.path.join(self.config.PACKAGE_DIR, 'whole_transcriptome/expression/deseq2batch.r')
        }
        self.file = {
            'json': os.path.join(self.work_dir, 'input.json'),
            'count': os.path.join(self.work_dir, 'count.txt'),
            'design': os.path.join(self.work_dir, 'design.txt'),
            'group': os.path.join(self.work_dir, 'group.txt'),
            'contrast': os.path.join(self.work_dir, 'contrast.txt'),
            'pheno': os.path.join(self.work_dir, 'pheno.txt')
        }

    def run(self):
        super(Deseq2batchTool, self).run()
        self.pre_deseq2batch()
        self.run_deseq2batch()
        self.set_output()
        self.end()

    def pre_deseq2batch(self):
        shutil.copy(self.option('group_table').path, self.file['group'])
        shutil.copy(self.option('control_table').path, self.file['contrast'])
        shutil.copy(self.option('pheno_matrix').path, self.file['pheno'])
        samples = list()
        groups = list()
        batchs = list()
        for line in open(self.file['group']):
            if line.strip() and line[0] != '#':
                sample, group = line.strip().split('\t')
                if sample not in samples:
                    samples.append(sample)
                    groups.append(group)
        for line in open(self.file['pheno']):
            if line.strip() and line[0] != '#':
                sample, batch = line.strip().split('\t')
                batchs.append(batch)

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
            'count': self.file['count'],
            'samples': samples,
            'groups': groups,
            'batchs': batchs,
            'contrasts': contrasts,
            'together': self.option('together'),
            'output': self.output_dir
        }, open(self.file['json'], 'w'), indent=4)

    def run_deseq2batch(self):
        cmd = '{} {} -i {}'.format(self.program['rscript'], self.script['deseq2batch'], self.file['json'])
        runcmd(self, 'run_deseq2batch', cmd)

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
            'id': 'deseq2batch_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'tool',
            'name': 'whole_transcriptome.expression.deseq2batch',
            'instant': False,
            'options': {
                'count_matrix': '/mnt/ilustre/users/sanger-dev/workspace/20191114/DiffExp_tsg_36088_6566_5962/ExpBuild/output/G.mRNA.all.txt',
                'group_table': '/mnt/ilustre/users/sanger-dev/workspace/20191114/DiffExp_tsg_36088_6566_5962/ExpBuild/output/group.txt',
                'control_table': '/mnt/ilustre/users/sanger-dev/workspace/20191114/DiffExp_tsg_36088_6566_5962/ExpBuild/output/control.txt',
                'pheno_matrix':'/mnt/ilustre/users/sanger-dev/workspace/20191114/DiffExp_tsg_36088_6566_5962/ExpBuild/output/pheno.txt',
                'together': 'True',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)

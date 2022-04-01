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


class CombatseqAgent(Agent):
    '''
    last_modify: 2019.12.12
    '''

    def __init__(self, parent):
        super(CombatseqAgent, self).__init__(parent)
        options = [
            {'name': 'count_matrix', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'group_table', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'batch_matrix', 'type': 'infile', 'format': 'whole_transcriptome.common'}


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
        super(CombatseqAgent, self).end()


class CombatseqTool(Tool):
    def __init__(self, config):
        super(CombatseqTool, self).__init__(config)
        self.program = {
            'python': 'program/Python/bin/python',
            'rscript': 'program/R-3.3.1/bin/Rscript',
        }
        self.script = {
            'combatseq': os.path.join(self.config.PACKAGE_DIR, 'tool_lab/Additive_model/batch/combatseq.r'),
        }
        self.file = {
            'json': os.path.join(self.work_dir, 'input.json'),
            'count': os.path.join(self.work_dir, 'count.txt'),
            'group': os.path.join(self.work_dir, 'group.txt'),
            'batch': os.path.join(self.work_dir, 'batch.txt'),
        }

    def run(self):
        super(CombatseqTool, self).run()
        self.pre()
        self.run_combatseq()
        self.set_output()
        self.end()

    def pre(self):
        shutil.copy(self.option('count_matrix').path, self.file['count'])
        shutil.copy(self.option('group_table').path, self.file['group'])
        shutil.copy(self.option('batch_matrix').path, self.file['batch'])
        samples = list()
        groups = list()
        batchs = list()
        for line in open(self.file['group']):
            if line.strip() and line[0] != '#':
                sample, group = line.strip().split('\t')
                if sample not in samples:
                    samples.append(sample)
                    groups.append(group)
        for line in open(self.file['batch']):
            if line.strip() and line[0] != '#':
                sample, batch = line.strip().split('\t')
                batchs.append(batch)
        df = pd.read_table(self.option('count_matrix').path, index_col=0)
        df = df.reindex(samples, axis=1)
        df = df.astype(int)
        df.to_csv(self.file['count'], sep='\t')
        json.dump({
            'output': self.output_dir,
            'count': self.file['count'],
            'groups': groups,
            'batchs':batchs
        }, open(self.file['json'], 'w'), indent=4)

    def run_combatseq(self):
        cmd = '{} {} -i {}'.format(self.program['rscript'], self.script['combatseq'], self.file['json'])
        runcmd(self, 'run_sva', cmd)



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
            'id': 'combatseq_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'tool',
            'name': 'tool_lab.Additive_model.batch.combatseq',
            'instant': False,
            'options': {
                'count_matrix': '/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/count_MJ2018.txt',
                'group_table': '/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/group.txt',
                'batch_matrix': '/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/tpm_batch.txt'
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)

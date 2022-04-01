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


class BatchAgent(Agent):
    '''
    last_modify: 2019.12.12
    '''

    def __init__(self, parent):
        super(BatchAgent, self).__init__(parent)
        options = [
            {'name': 'count_matrix', 'type': 'infile', 'format': 'ref_rna_v2.common'},
            {'name': 'group_table', 'type': 'infile', 'format': 'ref_rna_v2.common'},
            {'name': 'batch_matrix', 'type': 'infile', 'format': 'ref_rna_v2.common'},
            {'name': 'batch_method', 'type': 'string'},
            {'name': 'count_batch', 'type': 'outfile', 'format': 'ref_rna_v2.common'}

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
        super(BatchAgent, self).end()


class BatchTool(Tool):
    def __init__(self, config):
        super(BatchTool, self).__init__(config)
        software_dir = self.config.SOFTWARE_DIR
        self.gcc = software_dir + '/gcc/5.1.0/bin'
        self.gcc_lib = software_dir + '/gcc/5.1.0/lib64'
        self.set_environ(PATH=self.gcc, LD_LIBRARY_PATH=self.gcc_lib)
        self.program = {
            'python': 'program/Python/bin/python',
            'rscript': 'program/R-3.3.1/bin/Rscript',
        }
        self.script = {
            'combat': os.path.join(self.config.PACKAGE_DIR, 'medical_transcriptome/batch/combat.r'),
            'removebtacheffect': os.path.join(self.config.PACKAGE_DIR, 'medical_transcriptome/batch/removebatcheffect.r')
        }
        self.file = {
            'json': os.path.join(self.work_dir, 'input.json'),
            'count': os.path.join(self.work_dir, 'count.txt'),
            'group': os.path.join(self.work_dir, 'group.txt'),
            'batch': os.path.join(self.work_dir, 'batch.txt'),
            'count_batch': os.path.join(self.output_dir, 'count_batch.txt')
        }

    def run(self):
        super(BatchTool, self).run()
        self.pre()
        if self.option('batch_method') == 'combat':
            self.run_combat()
        if self.option('batch_method') == 'removeBatchEffect':
            self.run_removebatcheffect()
        self.set_output()
        self.end()

    def pre(self):
        shutil.copy(self.option('count_matrix').path, self.file['count'])
        shutil.copy(self.option('group_table').path, self.file['group'])
        shutil.copy(self.option('batch_matrix').path, self.file['batch'])
        batchs = list()
        samples = list()
        groups = list()
        sample_batch = dict()
        for line in open(self.file['group']):
            if line.strip() and line[0] != '#':
                sample, group = line.strip().split('\t')
                if sample not in samples:
                    samples.append(sample)
                    groups.append(group)

        for line in open(self.file['batch']):
            if line.strip() and line[0] != '#':
                sample, batch = line.strip().split('\t')
                sample_batch[sample] = batch
        for s in samples:
            batchs.append(sample_batch[s])
        df = pd.read_table(self.option('count_matrix').path, index_col=0)
        df = df.reindex(samples, axis=1)
        # df = df.astype(int)
        df.to_csv(self.file['count'], sep='\t')
        json.dump({
            'output': self.output_dir,
            'count': self.file['count'],
            'groups': groups,
            'batchs': batchs,
        }, open(self.file['json'], 'w'), indent=4)

    def run_combat(self):
        cmd = '{} {} -i {}'.format(self.program['rscript'], self.script['combat'], self.file['json'])
        runcmd(self, 'run_combat', cmd)

    def run_removebatcheffect(self):
        cmd = '{} {} -i {}'.format(self.program['rscript'], self.script['removebtacheffect'], self.file['json'])
        runcmd(self, 'run_removebatcheffect', cmd)


    def set_output(self):
        self.option('count_batch').set_path(self.file['count_batch'])



class TestFunction(unittest.TestCase):
    '''
    This is test for the tool. Just run this script to do test.
    '''

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'batch_limma{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'tool',
            'name': 'ref_rna_v3.batch.batch',
            'instant': False,
            'options': {
                'count_matrix': '/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/ref_rna_v3/batch/unigene.tpm.matrix.annot.txt',
                'group_table': '/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/ref_rna_v3/batch/group_table.txt',
                'batch_matrix': '/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/ref_rna_v3/batch/batch_new.txt',
                'batch_method': 'combat'
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)

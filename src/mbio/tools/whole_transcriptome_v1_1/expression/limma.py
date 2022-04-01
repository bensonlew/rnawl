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


class LimmaAgent(Agent):
    '''
    last_modify: 2019.11.15
    '''

    def __init__(self, parent):
        super(LimmaAgent, self).__init__(parent)
        options = [
            {'name': 'count_matrix', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'group_table', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'control_table', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'is_batch', 'type': 'bool', 'default': False},
            {'name': 'has_batch', 'type': 'bool', 'default': True},
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
        super(LimmaAgent, self).end()


class LimmaTool(Tool):
    def __init__(self, config):
        super(LimmaTool, self).__init__(config)
        software_dir = self.config.SOFTWARE_DIR
        self.gcc = software_dir + '/gcc/5.1.0/bin'
        self.gcc_lib = software_dir + '/gcc/5.1.0/lib64'
        self.set_environ(PATH=self.gcc, LD_LIBRARY_PATH=self.gcc_lib)
        self.r_path = software_dir + "/program/R-3.3.1/bin:$PATH"
        self._r_home = software_dir + "/program/R-3.3.1/lib64/R/"
        self._LD_LIBRARY_PATH = software_dir + "/program/R-3.3.1/lib64/R/lib:$LD_LIBRARY_PATH"
        self.set_environ(PATH=self.r_path, R_HOME=self._r_home, LD_LIBRARY_PATH=self._LD_LIBRARY_PATH)
        self.program = {
            'rscript': '/program/R-3.3.1/bin/Rscript',
        }
        self.script = {
            'limma': os.path.join(self.config.PACKAGE_DIR, 'whole_transcriptome_v1_1/expression/limma.r')
        }
        self.file = {
            'json': os.path.join(self.work_dir, 'input.json'),
            'count': os.path.join(self.work_dir, 'count.txt'),
            'group': os.path.join(self.work_dir, 'group.txt'),
            'control': os.path.join(self.work_dir, 'control.txt')
        }

    def run(self):
        super(LimmaTool, self).run()
        self.pre_limma()
        self.run_limma()
        self.set_output()
        self.end()

    def pre_limma(self):
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
        if self.option('is_batch') == True:
            batchs_dict = dict()
            for line in open(self.option('batch_matrix').path):
                if line.strip() and line[0] != '#':
                    sample, batch = line.strip().split('\t')
                    batchs_dict[sample] = batch
            batchs = [batchs_dict[i] for i in samples]
        contrasts = list()
        for line in open(self.file['control']):
            if line.strip() and line[0] != '#':
                ctrl, case = line.strip().split('\t')
                contrasts.append({'ctrl': ctrl, 'case': case})
        df = pd.read_table(self.option('count_matrix').path, index_col=0)
        df = df.reindex(samples, axis=1)
        df.to_csv(self.file['count'], sep='\t')
        if len(set(samples)) == len(set(groups)):
            singal_sample = 'yes'
        else:
            singal_sample = 'no'
        opts = {
            'count': self.file['count'],
            'groups': groups,
            'contrasts': contrasts,
            'output': self.output_dir,
            'is_batch': self.option('is_batch'),
            'singal_sample': singal_sample,
        }
        if self.option('is_batch') == True:
            opts.update({
                'batchs': batchs
            })
        json.dump(opts, open(self.file['json'], 'w'), indent=4)

    def run_limma(self):
        cmd = '{} {} -i {}'.format(self.program['rscript'], self.script['limma'], self.file['json'])
        runcmd(self, 'run_limma', cmd)

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
            'id': 'limma_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'tool',
            'name': 'whole_transcriptome.expression.limma',
            'instant': False,
            'options': {
                'count_matrix': '/mnt/ilustre/users/sanger-dev/workspace/20201028/DiffExp_tsg_219018_4794_99/DiffExp/Edger/count.txt',
                'group_table': '/mnt/ilustre/users/sanger-dev/workspace/20201028/DiffExp_tsg_219018_4794_99/DiffExp/Edger/group.txt',
                'control_table': '/mnt/ilustre/users/sanger-dev/workspace/20201028/DiffExp_tsg_219018_4794_99/DiffExp/Edger/control.txt'
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)

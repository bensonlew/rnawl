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
            {'name': 'exp_matrix', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'pheno_matrix', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'batch', 'type': 'string'}

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
        self.program = {
            'python': 'program/Python/bin/python',
            'rscript': 'program/R-3.3.1/bin/Rscript',
        }
        self.script = {
            'combat': os.path.join(self.config.PACKAGE_DIR, 'whole_transcriptome/expression/combat.r'),
            'removebtacheffect': os.path.join(self.config.PACKAGE_DIR, 'whole_transcriptome/expression/limma.r')
        }
        self.file = {
            'json': os.path.join(self.work_dir, 'input.json'),
            'exp': os.path.join(self.work_dir, 'exp_matrix.txt'),
            'pheno': os.path.join(self.work_dir, 'pheno_matrix.txt')
        }

    def run(self):
        super(BatchTool, self).run()
        self.pre()
        if self.option('batch') == 'combat':
            self.run_combat()
        if self.option('batch') == 'limma':
            self.run_removebatcheffect()
        self.set_output()
        self.end()

    def pre(self):
        shutil.copy(self.option('exp_matrix').path, self.file['exp'])
        shutil.copy(self.option('pheno_matrix').path, self.file['pheno'])
        json.dump({
            'output': self.output_dir,
            'exp': self.file['exp'],
            'pheno': self.file['pheno']
        }, open(self.file['json'], 'w'), indent=4)

    def run_combat(self):
        cmd = '{} {} -i {}'.format(self.program['rscript'], self.script['combat'], self.file['json'])
        runcmd(self, 'run_combat', cmd)

    def run_removebatcheffect(self):
        cmd = '{} {} -i {}'.format(self.program['rscript'], self.script['removebtacheffect'], self.file['json'])
        runcmd(self, 'run_removebatcheffect', cmd)


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
            'id': 'batch_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'tool',
            'name': 'whole_transcriptome.expression.batch',
            'instant': False,
            'options': {
                'exp_matrix': '/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/batch/r_test.txt',
                'pheno_matrix': '/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/batch/pheno_test.txt',
                'batch': 'limma'
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)

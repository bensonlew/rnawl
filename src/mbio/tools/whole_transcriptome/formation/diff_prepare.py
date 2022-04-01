# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import json
import os
import unittest

from biocluster.agent import Agent
from biocluster.tool import Tool

from mbio.packages.whole_transcriptome.utils import runcmd


class DiffPrepareAgent(Agent):
    '''
    last_modify: 2019.11.14
    '''

    def __init__(self, parent):
        super(DiffPrepareAgent, self).__init__(parent)
        LEVEL = ('G', 'T')
        CATEGORY = ('mRNA', 'lncRNA', 'miRNA', 'circRNA')
        KIND = ('all', 'ref')
        options = [
            {'name': 'task_id', 'type': 'string', 'default': None},
            {'name': 'level', 'type': 'string', 'default': LEVEL[1]},
            {'name': 'category', 'type': 'string', 'default': CATEGORY[0]},
            {'name': 'kind', 'type': 'string', 'default': KIND[0]},
            {'name': 'background', 'type': 'string', 'default': None},
            {'name': 'count_matrix', 'type': 'outfile', 'format': 'whole_transcriptome.common'},
            {'name': 'kind_table', 'type': 'outfile', 'format': 'whole_transcriptome.common'},
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
        super(DiffPrepareAgent, self).end()


class DiffPrepareTool(Tool):
    def __init__(self, config):
        super(DiffPrepareTool, self).__init__(config)
        self.program = {
            'python': 'program/Python/bin/python'
        }
        self.script = {
            'diff_prepare': os.path.join(self.config.PACKAGE_DIR, 'whole_transcriptome/formation/diff_prepare.py')
        }
        args = [self.option('level'), self.option('category'), self.option('kind')]
        self.file = {
            'json': os.path.join(self.work_dir, 'input.json'),
            'count_matrix': os.path.join(self.output_dir, '{}.{}.{}.txt'.format(*args)),
            'kind_table': os.path.join(self.output_dir, 'kind.txt')
        }

    def run(self):
        super(DiffPrepareTool, self).run()
        self.pre_diff_prepare()
        self.run_diff_prepare()
        self.set_output()
        self.end()

    def pre_diff_prepare(self):
        json.dump({
            'task_id': self.option('task_id'),
            'level': self.option('level'),
            'category': self.option('category'),
            'kind': self.option('kind'),
            'background': self.option('background'),
            'work_dir': self.work_dir,
            'count_matrix': self.file['count_matrix'],
            'kind_table': self.file['kind_table']
        }, open(self.file['json'], 'w'), indent=4)

    def run_diff_prepare(self):
        cmd = '{} {} -j {}'.format(self.program['python'], self.script['diff_prepare'], self.file['json'])
        if self.config.DBVersion == 1:
            DBVersion =1
        else:
            DBVersion = 0
        cmd += ' -v {}'.format(DBVersion)
        runcmd(self, 'run_diff_prepare', cmd)

    def set_output(self):
        self.option('count_matrix').set_path(self.file['count_matrix'])
        self.option('kind_table').set_path(self.file['kind_table'])


class TestFunction(unittest.TestCase):
    '''
    This is test for the tool. Just run this script to do test.
    '''

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'diff_prepare_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'tool',
            'name': 'whole_transcriptome.formation.diff_prepare',
            'instant': False,
            'options': {
                'task_id': 'tsg_36088',
                'level': 'T',
                'category': 'mRNA',
                'kind': 'ref',
                'background': 'mRNA,lncRNA'
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)

# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import json
import os
import unittest

from biocluster.agent import Agent
from biocluster.tool import Tool

from mbio.packages.whole_transcriptome.utils import runcmd


class ExpBuildAgent(Agent):
    '''
    last_modify: 2019.11.13
    '''

    def __init__(self, parent):
        super(ExpBuildAgent, self).__init__(parent)
        LIBRARY = ('long', 'small', 'circle')
        LEVEL = ('G', 'T')
        CATEGORY = ('mRNA', 'lncRNA', 'miRNA', 'circRNA')
        KIND = ('all', 'ref', 'new')
        options = [
            {'name': 'task_id', 'type': 'string', 'default': None},
            {'name': 'library', 'type': 'string', 'default': None},
            {'name': 'level', 'type': 'string', 'default': LEVEL[1]},
            {'name': 'category', 'type': 'string', 'default': CATEGORY[0]},
            {'name': 'kind', 'type': 'string', 'default': KIND[0]},
            {'name': 'group_dict', 'type': 'string', 'default': None},
            {'name': 'control_id', 'type': 'string', 'default': None},
            {'name': 'exp_matrix', 'type': 'outfile', 'format': 'whole_transcriptome.common'},
            {'name': 'group_table', 'type': 'outfile', 'format': 'whole_transcriptome.common'},
            {'name': 'control_table', 'type': 'outfile', 'format': 'whole_transcriptome.common'},
            {'name': 'is_rmbe', 'type': 'string', 'default': 'false'}
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
        super(ExpBuildAgent, self).end()


class ExpBuildTool(Tool):
    def __init__(self, config):
        super(ExpBuildTool, self).__init__(config)
        self.program = {
            'python': 'program/Python/bin/python'
        }
        self.script = {
            'exp_build': os.path.join(self.config.PACKAGE_DIR, 'whole_transcriptome_v1_1/formation/exp_build.py')
        }
        if self.option('library'):
            args = [self.option('library'), self.option('level'), self.option('kind')]
        else:
            args = [self.option('level'), self.option('category'), self.option('kind')]
        self.file = {
            'json': os.path.join(self.work_dir, 'input.json'),
            'exp_matrix': os.path.join(self.output_dir, '{}.{}.{}.txt'.format(*args)),
            'group_table': os.path.join(self.output_dir, 'group.txt'),
            'control_table': os.path.join(self.output_dir, 'control.txt')
        }

    def run(self):
        super(ExpBuildTool, self).run()
        if self.option('library'):
            self.pre_exp_build_library()
        else:
            self.pre_exp_build_category()
        self.run_exp_build()
        self.set_output()
        self.end()

    def pre_exp_build_library(self):
        json.dump({
            'task_id': self.option('task_id'),
            'library': self.option('library'),
            'level': self.option('level'),
            'kind': self.option('kind'),
            'group_dict': self.option('group_dict'),
            'control_id': self.option('control_id'),
            'exp_matrix': self.file['exp_matrix'],
            'group_table': self.file['group_table'],
            'control_table': self.file['control_table'],
            'is_rmbe': self.option('is_rmbe')
        }, open(self.file['json'], 'w'), indent=4)

    def pre_exp_build_category(self):
        json.dump({
            'task_id': self.option('task_id'),
            'level': self.option('level'),
            'category': self.option('category'),
            'kind': self.option('kind'),
            'group_dict': self.option('group_dict'),
            'control_id': self.option('control_id'),
            'exp_matrix': self.file['exp_matrix'],
            'group_table': self.file['group_table'],
            'control_table': self.file['control_table'],
            'is_rmbe': self.option('is_rmbe')
        }, open(self.file['json'], 'w'), indent=4)

    def run_exp_build(self):
        cmd = '{} {} -j {}'.format(self.program['python'], self.script['exp_build'], self.file['json'])
        if self.config.DBVersion == 1:
            DBVersion = 1
        else:
            DBVersion = 0
        cmd += ' -v {}'.format(DBVersion)
        runcmd(self, 'run_exp_build', cmd)
        # cmd = '{} {} {}'.format(self.program['python'], self.script['exp_build'], self.file['json'])
        # runcmd(self, 'run_exp_build', cmd)

    def set_output(self):
        self.option('exp_matrix').set_path(self.file['exp_matrix'])
        self.option('group_table').set_path(self.file['group_table'])
        if self.option('control_id'):
            self.option('control_table').set_path(self.file['control_table'])


class TestFunction(unittest.TestCase):
    '''
    This is test for the tool. Just run this script to do test.
    '''

    def test_library(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'exp_build_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'tool',
            'name': 'whole_transcriptome.formation.exp_build',
            'instant': False,
            'options': {
                'task_id': 'tsg_36088',
                'library': 'long',
                'level': 'T',
                'kind': 'ref',
                'group_dict': json.dumps({'Acute': ['A1', 'A2'], 'Stable': ['S1', 'S3'], 'control': ['C2', 'C3']}),
                'control_id': '5dcb545a17b2bf08b8f25cee'
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

    def test_category(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'exp_build_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'tool',
            'name': 'whole_transcriptome.formation.exp_build',
            'instant': False,
            'options': {
                'task_id': 'tsg_36088',
                'level': 'T',
                'category': 'mRNA',
                'kind': 'ref',
                'group_dict': json.dumps({'Acute': ['A1', 'A2'], 'Stable': ['S1', 'S3'], 'control': ['C2', 'C3']}),
                'control_id': '5dcb545a17b2bf08b8f25cee'
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test_category')])
    unittest.TextTestRunner(verbosity=2).run(suite)

# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import os
import shutil
import unittest

from biocluster.agent import Agent
from biocluster.tool import Tool

from mbio.packages.denovo_rna_v2.functions import reminder, runcmd


class AlignBeginAgent(Agent):
    '''
    last_modify: 2019.08.07
    '''

    def __init__(self, parent):
        super(AlignBeginAgent, self).__init__(parent)
        PURPOSE_OPTIONS = ('query', 'subject')
        UNIT_OPTIONS = ('nucl', 'prot')
        EXP_LEVEL_OPTIONS = ('T', 'G')
        SOURCE_OPTIONS = ('filter', 'raw')
        options = [
            {'name': 'purpose', 'type': 'string', 'default': PURPOSE_OPTIONS[0]},
            {'name': 'unit_type', 'type': 'string', 'default': UNIT_OPTIONS[0]},
            {'name': 'task_id', 'type': 'string', 'default': None},
            {'name': 'exp_level', 'type': 'string', 'default': EXP_LEVEL_OPTIONS[0]},
            {'name': 'geneset_id', 'type': 'string', 'default': 'All'},
            {'name': 'source', 'type': 'string', 'default': SOURCE_OPTIONS[0]},
            {'name': 'denovo', 'type': 'outfile', 'format': 'denovo_rna_v2.fasta'},
        ]
        self.add_option(options)

    @reminder
    def check_options(self):
        for k, v in self._options.items():
            self.logger.debug('{} = {}'.format(k, v.value))

    @reminder
    def set_resource(self):
        self._cpu = 1
        self._memory = '4G'

    @reminder
    def end(self):
        super(AlignBeginAgent, self).end()


class AlignBeginTool(Tool):
    def __init__(self, config):
        super(AlignBeginTool, self).__init__(config)
        self.program = {
            'python': 'miniconda2/bin/python'
        }
        self.script = {
            'get_denovo_seq': os.path.join(self.config.PACKAGE_DIR, 'denovo_rna_v2/alignment/get_denovo_seq.py')
        }
        self.file = {
            'query': os.path.join(self.output_dir, 'query.{}.{}.{}.fasta'.format(self.option('unit_type'), self.option('exp_level'), self.option('geneset_id'))),
            'subject': os.path.join(self.output_dir, 'subject.{}.fasta'.format(self.option('source'))),
        }

    @reminder
    def run(self):
        super(AlignBeginTool, self).run()
        self.run_get_denovo_seq()
        self.set_output()
        self.end()

    @reminder
    def run_get_denovo_seq(self):
        cmd = '{} {} {}'.format(self.program['python'], self.script['get_denovo_seq'], self.option('purpose'))
        cmd += ' --task {}'.format(self.option('task_id'))
        if self.option('purpose') == 'query':
            cmd += ' --unit {}'.format(self.option('unit_type'))
            cmd += ' --level {}'.format(self.option('exp_level'))
            cmd += ' --geneset {}'.format(self.option('geneset_id'))
        elif self.option('purpose') == 'subject':
            cmd += ' --source {}'.format(self.option('source'))
        if self.config.DBVersion == 1:
            DBVersion = 1
        else:
            DBVersion = 0
        cmd += ' -v {}'.format(DBVersion)
        cmd += ' --output {}'.format(self.file[self.option('purpose')])
        cmd_name = 'run_get_denovo_seq'
        runcmd(self, cmd_name, cmd)

    @reminder
    def set_output(self):
        self.option('denovo').set_path(self.file[self.option('purpose')])


class TestFunction(unittest.TestCase):
    '''
    This is test for the tool. Just run this script to do test.
    '''

    def test_query_nucl_G(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'align_begin_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'tool',
            'name': 'denovo_rna_v2.alignment.align_begin',
            'instant': False,
            'options': {
                'purpose': 'query',
                'unit_type': 'nucl',
                'task_id': 'tsg_33857',
                'exp_level': 'G',
                'geneset_id': 'All',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

    def test_query_nucl_T(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'align_begin_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'tool',
            'name': 'denovo_rna_v2.alignment.align_begin',
            'instant': False,
            'options': {
                'purpose': 'query',
                'unit_type': 'nucl',
                'task_id': 'tsg_33857',
                'exp_level': 'T',
                'geneset_id': '5cb64c0d17b2bf6187b11fbd',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

    def test_query_prot_G(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'align_begin_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'tool',
            'name': 'denovo_rna_v2.alignment.align_begin',
            'instant': False,
            'options': {
                'purpose': 'query',
                'unit_type': 'prot',
                'task_id': 'tsg_33857',
                'exp_level': 'G',
                'geneset_id': '5d47b73217b2bf35cade6f7b',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

    def test_subject_filter(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'align_begin_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'tool',
            'name': 'denovo_rna_v2.alignment.align_begin',
            'instant': False,
            'options': {
                'purpose': 'subject',
                'task_id': 'tsg_33857',
                'source': 'filter'
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test_subject_filter')])
    unittest.TextTestRunner(verbosity=2).run(suite)

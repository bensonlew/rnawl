# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import os
import unittest

from biocluster.agent import Agent
from biocluster.tool import Tool

from mbio.packages.whole_transcriptome.utils import runcmd


class DiffSplitAgent(Agent):
    '''
    last_modify: 2019.11.08
    '''

    def __init__(self, parent):
        super(DiffSplitAgent, self).__init__(parent)
        options = [
            {'name': 'diff_dir', 'type': 'infile', 'format': 'whole_transcriptome.common_dir'},
            {'name': 'relation', 'type': 'infile', 'format': 'whole_transcriptome.common'},
        ]
        self.add_option(options)

    def check_options(self):
        for k, v in self._options.items():
            self.logger.debug('{} -> {}'.format(k, v.value))
        else:
            return True

    def set_resource(self):
        self._cpu = 1
        self._memory = '5G'

    def end(self):
        super(DiffSplitAgent, self).end()


class DiffSplitTool(Tool):
    def __init__(self, config):
        super(DiffSplitTool, self).__init__(config)
        self.program = {
            'python': 'miniconda2/bin/python'
        }
        self.script = {
            'diff_split': os.path.join(self.config.PACKAGE_DIR, 'whole_transcriptome/formation/diff_split.py')
        }

    def run(self):
        super(DiffSplitTool, self).run()
        self.run_diff_split()
        self.set_output()
        self.end()

    def run_diff_split(self):
        cmd = '{} {}'.format(self.program['python'], self.script['diff_split'])
        cmd += ' -i {}'.format(self.option('diff_dir').path)
        cmd += ' -t {}'.format(self.option('relation').path)
        cmd += ' -o {}'.format(self.output_dir)
        runcmd(self, 'run_diff_split', cmd)

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
            'id': 'diff_split_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'tool',
            'name': 'whole_transcriptome.formation.diff_split',
            'instant': False,
            'options': {
                'diff_dir': '/mnt/ilustre/users/sanger-dev/workspace/20191106/WholeTranscriptome_tsg_36088/Transfer/output/diff_exp_t',
                'relation': '/mnt/ilustre/users/sanger-dev/workspace/20191106/WholeTranscriptome_tsg_36088/Transfer/output/large_gush/filter_by_express/filtered_file/trans_type.xls'
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)

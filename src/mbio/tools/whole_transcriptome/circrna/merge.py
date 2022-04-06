# -*- coding: utf-8 -*-
# __author__ = 'zoujiaxun'

import os
import unittest

from biocluster.agent import Agent
from biocluster.tool import Tool

from mbio.packages.whole_transcriptome.utils import runcmd


class MergeAgent(Agent):
    '''
    last_modify: 2019.10.22
    '''

    def __init__(self, parent):
        super(MergeAgent, self).__init__(parent)
        options = [
            {'name': 'sample', 'type': 'string', 'default': None},
            {'name': 'find_circ_filter', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'ciri_circ_filter', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'circRNA_merge', 'type': 'outfile', 'format': 'whole_transcriptome.common'}
        ]
        self.add_option(options)

    def check_options(self):
        pass

    def set_resource(self):
        self._cpu = 1
        self._memory = '4G'

    def end(self):
        super(MergeAgent, self).end()


class MergeTool(Tool):
    def __init__(self, config):
        super(MergeTool, self).__init__(config)
        self.program = {
            # 'python': 'miniconda2/bin/python',
            'python': 'bioinfo/rna/miniconda2/bin/python'
        }
        self.script = {
            'merge': os.path.join(self.config.PACKAGE_DIR, 'whole_transcriptome/circrna/merge.py')
        }
        self.file = {
            'circRNA_merge': os.path.join(self.output_dir, '{}_circRNA_merge.txt'.format(self.option('sample'))),
        }

    def run(self):
        super(MergeTool, self).run()
        self.run_merge()
        self.set_output()
        self.end()

    def run_merge(self):
        cmd = '{} {}'.format(self.program['python'], self.script['merge'])
        cmd += ' -f {} -c {} -o {}'.format(self.option('find_circ_filter').path, self.option('ciri_circ_filter').path,
                                           self.file['circRNA_merge'])
        cmd_name = 'run_merge'
        runcmd(self, cmd_name, cmd)

    def set_output(self):
        self.option('circRNA_merge').set_path(self.file['circRNA_merge'])


class TestFunction(unittest.TestCase):
    '''
    This is test for the tool. Just run this script to do test.
    '''

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'merge_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'tool',
            'name': 'whole_transcriptome.circrna.merge',
            'instant': False,
            'options': {
                'find_circ_filter': '/mnt/ilustre/users/sanger-dev/workspace/20191021/Single_circ_brush_7068_2074/CircBrush/Ciriaddfindcirc1/Findcircbwa/output/S1_find_circ_filter.txt',
                'ciri_circ_filter': '/mnt/ilustre/users/sanger-dev/workspace/20191021/Single_circ_brush_7068_2074/CircBrush/Ciriaddfindcirc1/Ciri2/output/S1_ciri_circ_filter.txt',
                'sample': 'S1'
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)

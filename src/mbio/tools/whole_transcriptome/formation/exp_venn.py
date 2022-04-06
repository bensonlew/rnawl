# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import os
import unittest

from biocluster.agent import Agent
from biocluster.tool import Tool

from mbio.packages.whole_transcriptome.utils import runcmd


class ExpVennAgent(Agent):
    '''
    last_modify: 2019.09.20
    '''

    def __init__(self, parent):
        super(ExpVennAgent, self).__init__(parent)
        options = [
            {'name': 'exp_matrix', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'group_table', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'threshold', 'type': 'float', 'default': 1.0},
            {'name': 'venn_table', 'type': 'outfile', 'format': 'whole_transcriptome.common'},
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
        super(ExpVennAgent, self).end()


class ExpVennTool(Tool):
    def __init__(self, config):
        super(ExpVennTool, self).__init__(config)
        self.program = {
            'python': 'miniconda2/bin/python'
        }
        self.script = {
            'exp_venn': os.path.join(self.config.PACKAGE_DIR, 'whole_transcriptome/formation/exp_venn.py')
        }
        self.file = {
            'venn': os.path.join(self.output_dir, 'venn.txt')
        }

    def run(self):
        super(ExpVennTool, self).run()
        self.run_exp_venn()
        self.set_output()
        self.end()

    def run_exp_venn(self):
        cmd = '{} {}'.format(self.program['python'], self.script['exp_venn'])
        cmd += ' -i {}'.format(self.option('exp_matrix').path)
        cmd += ' -g {}'.format(self.option('group_table').path)
        cmd += ' -t {}'.format(self.option('threshold'))
        cmd += ' -o {}'.format(self.file['venn'])
        runcmd(self, 'run_exp_venn', cmd)

    def set_output(self):
        self.option('venn_table').set_path(self.file['venn'])


class TestFunction(unittest.TestCase):
    '''
    This is test for the tool. Just run this script to do test.
    '''

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'exp_venn_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'tool',
            'name': 'whole_transcriptome.formation.exp_venn',
            'instant': False,
            'options': {
                # 'exp_matrix': '/mnt/ilustre/users/sanger-dev/workspace/20190919/Single_exp_make_9498_2417/ExpMake/output/mrna/G.tpm.txt',
                'exp_matrix': '/mnt/ilustre/users/sanger-dev/workspace/20190919/Single_exp_make_9498_2417/ExpMake/output/mrna/T.tpm.txt',
                'group_table': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/whole_transcriptome/MJ20190118060/largeRNA/example_group.txt',
                'threshold': 2.0,
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)

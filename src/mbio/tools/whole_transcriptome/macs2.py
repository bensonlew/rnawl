# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import os
import unittest

from biocluster.agent import Agent
from biocluster.tool import Tool

from mbio.packages.whole_transcriptome.utils import runcmd


class Macs2Agent(Agent):
    '''
    last_modify: 2020.02.17
    '''

    def __init__(self, parent):
        super(Macs2Agent, self).__init__(parent)
        options = [
            {'name': 'treatment', 'type': 'infile', 'format': 'whole_transcriptome.bam'},
            {'name': 'control', 'type': 'infile', 'format': 'whole_transcriptome.bam'},
            {'name': 'format', 'type': 'string', 'default': 'BAMPE'},
            {'name': 'gsize', 'type': 'string', 'default': 'hs'},
            {'name': 'name', 'type': 'string', 'default': 'NA'},
            {'name': 'bdg', 'type': 'bool', 'default': True},
            {'name': 'verbose', 'type': 'int', 'default': 2},
            {'name': 'qvalue', 'type': 'float', 'default': 0.05},
        ]
        self.add_option(options)

    def check_options(self):
        for k, v in self._options.items():
            self.logger.debug('{} => {}'.format(k, v.value))

    def set_resource(self):
        self._cpu = 1
        self._memory = '20G'

    def end(self):
        super(Macs2Agent, self).end()


class Macs2Tool(Tool):
    def __init__(self, config):
        super(Macs2Tool, self).__init__(config)
        self.program = {
            'macs2': 'bioinfo/medip/miniconda3/bin/macs2'
        }

    def run(self):
        super(Macs2Tool, self).run()
        self.run_macs2()
        self.set_output()
        self.end()

    def run_macs2(self):
        cmd = '{} callpeak'.format(self.program['macs2'])
        cmd += ' -t {}'.format(self.option('treatment').path)
        if self.option('control').is_set:
            cmd += ' -c {}'.format(self.option('control').path)
        cmd += ' -f {}'.format(self.option('format'))
        cmd += ' -g {}'.format(self.option('gsize'))
        cmd += ' --outdir {}'.format(self.output_dir)
        cmd += ' -n {}'.format(self.option('name'))
        cmd += ' --verbose {}'.format(self.option('verbose'))
        cmd += ' -q {}'.format(self.option('qvalue'))
        runcmd(self, 'run_macs2', cmd)

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
            'id': 'macs2_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'tool',
            'name': 'whole_transcriptome.macs2',
            'instant': False,
            'options': {
                'treatment': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/medip/mapping_data/YGXLTHZ1225-4-YG201902876-JJH_L2.bam',
                'name': 'YGXLTHZ1225-4-YG201902876-JJH_L2',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)

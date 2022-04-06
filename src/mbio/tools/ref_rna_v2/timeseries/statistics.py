# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

from biocluster.agent import Agent
from mbio.packages.ref_rna_v2.functions import toolfuncdeco, runcmd
from biocluster.tool import Tool
import os
import unittest

class StatisticsAgent(Agent):
    '''
    last_modify: 2019.06.18
    '''
    def __init__(self, parent):
        super(StatisticsAgent, self).__init__(parent)
        options = [
            {'name': 'method', 'type': 'string', 'default': 'SCM'},  # ['SCM', 'K']
            {'name': 'genetable', 'type': 'infile', 'format': 'ref_rna_v2.common'},
            {'name': 'profiletable', 'type': 'infile', 'format': 'ref_rna_v2.common'},
            {'name': 'kmeansclustertable', 'type': 'infile', 'format': 'ref_rna_v2.common'},
            {'name': 'matrix', 'type': 'infile', 'format': 'ref_rna_v2.matrix'},
            {'name': 'detailtable', 'type': 'outfile', 'format': 'ref_rna_v2.common'},
            {'name': 'clustertable', 'type': 'outfile', 'format': 'ref_rna_v2.common'},
        ]
        self.add_option(options)
        self.step.add_steps('statistics')
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.statistics.start()
        self.step.update()

    def stepfinish(self):
        self.step.statistics.finish()
        self.step.update()

    @toolfuncdeco
    def check_options(self):
        for k, v in self._options.items():
            self.logger.debug('{} = {}'.format(k, v.value))

    def set_resource(self):
        self._cpu = 1
        self._memory = '{}G'.format(os.path.getsize(self.option('genetable').path) / 1024 ** 3 + 8)

    @toolfuncdeco
    def end(self):
        super(StatisticsAgent, self).end()

class StatisticsTool(Tool):
    def __init__(self, config):
        super(StatisticsTool, self).__init__(config)
        self.python = 'miniconda2/bin/python'
        self.script = {
            'statistics': os.path.join(self.config.PACKAGE_DIR, 'ref_rna_v2/timeseries/statistics.py'),
        }
        self.file = {
            'detail': os.path.join(self.output_dir, 'detail.{}.tsv'.format(self.option('method'))),
            'cluster': os.path.join(self.output_dir, 'cluster.{}.tsv'.format(self.option('method')))
        }

    @toolfuncdeco
    def run(self):
        super(StatisticsTool, self).run()
        self.run_statistics()
        self.set_output()
        self.end()

    @toolfuncdeco
    def run_statistics(self):
        cmd = '{} {}'.format(self.python, self.script['statistics'])
        cmd += ' --method {}'.format(self.option('method'))
        cmd += ' --genetable {}'.format(self.option('genetable').path)
        if self.option('method') == 'SCM':
            cmd += ' -p {}'.format(self.option('profiletable').path)
        elif self.option('method') == 'K':
            cmd += ' -k {}'.format(self.option('kmeansclustertable').path)
        cmd += ' --matrix {}'.format(self.option('matrix').path)
        cmd += ' --output {}'.format(self.output_dir)
        cmd_name = 'run_statistics'
        runcmd(self, cmd_name, cmd)

    @toolfuncdeco
    def set_output(self):
        self.option('detailtable').set_path(self.file['detail'])
        self.option('clustertable').set_path(self.file['cluster'])

class TestFunction(unittest.TestCase):
    '''
    This is test for the tool. Just run this script to do test.
    '''
    def test_scm(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'statistics_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'tool',
            'name': 'ref_rna_v2.timeseries.statistics',
            'instant': False,
            'options': {
                'method': 'SCM',
                'genetable': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/ref_rna_v2/stem/output/setting_genetable.scm.txt',
                'profiletable': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/ref_rna_v2/stem/output/setting_profiletable.txt',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

    def test_k(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'statistics_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'tool',
            'name': 'ref_rna_v2.timeseries.statistics',
            'instant': False,
            'options': {
                'method': 'K',
                'genetable': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/ref_rna_v2/stem/output/setting_genetable.k.txt',
                'kmeansclustertable': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/ref_rna_v2/stem/output/setting_kmeansclustertable.txt',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test_scm'), TestFunction('test_k')])
    unittest.TextTestRunner(verbosity=2).run(suite)

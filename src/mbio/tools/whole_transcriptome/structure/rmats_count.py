# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

from biocluster.agent import Agent
from mbio.packages.whole_transcriptome.functions import toolfuncdeco, runcmd
from biocluster.tool import Tool
import os
import unittest

class RmatsCountAgent(Agent):
    '''
    last_modify: 2019.06.11
    '''
    def __init__(self, parent):
        super(RmatsCountAgent, self).__init__(parent)
        options = [
            {'name': 'jcpklist', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'jcecpklist', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'ojc', 'type': 'outfile', 'format': 'whole_transcriptome.common'},
            {'name': 'ojcec', 'type': 'outfile', 'format': 'whole_transcriptome.common'}
        ]
        self.add_option(options)
        self.step.add_steps('rmats_count')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    @toolfuncdeco
    def check_options(self):
        for k, v in self._options.items():
            self.logger.debug('{} = {}'.format(k, v.value))

    def step_start(self):
        self.step.rmats_count.start()
        self.step.update()

    def step_end(self):
        self.step.rmats_count.finish()
        self.step.update()

    def set_resource(self):
        self._cpu = 1
        self._memory = '10G'

    @toolfuncdeco
    def end(self):
        super(RmatsCountAgent, self).end()

class RmatsCountTool(Tool):
    def __init__(self, config):
        super(RmatsCountTool, self).__init__(config)
        self.program = {
            'python': 'program/Python/bin/python'
        }
        self.script = {
            'rmats_count': os.path.join(self.config.PACKAGE_DIR, 'ref_rna_v3/structure/rmats_count.py')
        }
        self.file = {
            'ojc': os.path.join(self.output_dir, 'sample.event.count.JC.txt'),
            'ojcec': os.path.join(self.output_dir, 'sample.event.count.JCEC.txt')
        }

    def run(self):
        super(RmatsCountTool, self).run()
        self.run_rmats_count()
        self.set_output()
        self.end()

    @toolfuncdeco
    def run_rmats_count(self):
        cmd = '{} {}'.format(self.program['python'], self.script['rmats_count'])
        cmd += ' --jc {}'.format(self.option('jcpklist').path)
        cmd += ' --jcec {}'.format(self.option('jcecpklist').path)
        cmd += ' --ojc {}'.format(self.file['ojc'])
        cmd += ' --ojcec {}'.format(self.file['ojcec'])
        cmd_name = 'run_rmats_count'
        runcmd(self, cmd_name, cmd)

    @toolfuncdeco
    def set_output(self):
        self.option('ojc').set_path(self.file['ojc'])
        self.option('ojcec').set_path(self.file['ojcec'])

class TestFunction(unittest.TestCase):
    '''
    This is test for the tool. Just run this script to do test.
    '''
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'rmats_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'tool',
            'name': 'ref_rna_v3.structure.rmats_count',
            'instant': False,
            'options': {
                'jcpklist': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/ref_rna_v2/rmats/count/jc.pk.list',
                'jcecpklist': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/ref_rna_v2/rmats/count/jcec.pk.list'
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)

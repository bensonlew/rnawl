# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

from biocluster.agent import Agent
from biocluster.tool import Tool
import os
import unittest

class FamilyPrepAgent(Agent):
    '''
    last_modify: 2019.04.19
    '''
    def __init__(self, parent):
        super(FamilyPrepAgent, self).__init__(parent)
        options = [
            {'name': 'lncrna_fa', 'type': 'infile', 'format': 'lnc_rna.fasta'},
            {'name': 'number', 'type': 'int', 'default': 10},
        ]
        self.add_option(options)
        self.step.add_steps('family_prep')
        self.on('start', self.step_start)
        self.on('end', self.step_finish)

    def step_start(self):
        self.step.family_prep.start()
        self.step.update()

    def step_finish(self):
        self.step.family_prep.finish()
        self.step.update()

    def check_options(self):
        self.logger.info('start check_options at {}'.format(self.__class__.__name__))
        for k, v in self._options.items():
            self.logger.debug('{} = {}'.format(k, v.value))
        self.logger.info('finish check_options at {}'.format(self.__class__.__name__))

    def set_resource(self):
        self._cpu = 1
        self._memory = '2G'

    def end(self):
        super(FamilyPrepAgent, self).end()

class FamilyPrepTool(Tool):
    def __init__(self, config):
        super(FamilyPrepTool, self).__init__(config)
        self.python = 'miniconda2/bin/python'
        self.family_prepare_py = os.path.join(self.config.PACKAGE_DIR, 'lnc_rna/family_prepare.py')

    def run(self):
        super(FamilyPrepTool, self).run()
        self.run_family_prep()
        self.set_output()
        self.end()

    def run_family_prep(self):
        cmd = '{} {}'.format(self.python, self.family_prepare_py)
        cmd += ' --input {}'.format(self.option('lncrna_fa').path)
        cmd += ' --number {}'.format(self.option('number'))
        cmd += ' --output {}'.format(self.output_dir)
        cmd_name = 'run_family_prep'
        self.run_code(cmd_name, cmd)

    def run_code(self, cmd_name, cmd, shell=False, block=True):
        if shell:
            cmd = self.config.SOFTWARE_DIR + '/' + cmd
        command = self.add_command(cmd_name, cmd, shell=shell)
        command.run()
        command.no_check = True
        if block:
            self.wait()
            for n, c in self.commands.items():
                if c.no_check:
                    if c.return_code == c.default_return_code:
                        c.no_check = False
                        self.logger.info('succeed in running {}'.format(n))
                    else:
                        self.set_error('fail to run {}, abord'.format(n))

    def set_output(self):
        self.logger.info('start set_output at {}'.format(self.__class__.__name__))
        self.logger.info('finish set_output at {}'.format(self.__class__.__name__))

class TestFunction(unittest.TestCase):
    '''
    This is test for the tool. Just run this script to do test.
    '''
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'family_prep_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'tool',
            'name': 'lnc_rna.structure.family_prep',
            'instant': False,
            'options': {
                'lncrna_fa': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/lnc_rna/lncrna_family/family_prep/all_lncrna.fa',
                'number': 10,
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)

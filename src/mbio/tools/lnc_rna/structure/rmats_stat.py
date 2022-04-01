# -*- coding: utf-8 -*-
# __author__ = 'jinlinfang, qinjincheng'

from biocluster.agent import Agent
from biocluster.tool import Tool
import os
import shutil
import unittest

class RmatsStatAgent(Agent):
    '''
    last_modify: 2019.03.20
    '''
    def __init__(self, parent):
        super(RmatsStatAgent, self).__init__(parent)
        options = [
            {'name': 'result_dir', 'type': 'infile', 'format': 'lnc_rna.common_dir'},
            {'name': 'delta_psi', 'type': 'float', 'default': None},
            {'name': 'significant_diff', 'type': 'string', 'default': None},
            {'name': 'significant_value', 'type': 'float', 'default': None},
        ]
        self.add_option(options)
        self.step.add_steps('rmats_stat')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def check_options(self):
        self.logger.info('start check_options at {}'.format(self.__class__.__name__))
        for k, v in self._options.items():
            self.logger.debug('{} = {}'.format(k, v.value))
        self.logger.info('finish check_options at {}'.format(self.__class__.__name__))

    def step_start(self):
        self.step.rmats_stat.start()
        self.step.update()

    def step_end(self):
        self.step.rmats_stat.finish()
        self.step.update()

    def set_resource(self):
        self._cpu = 1
        self._memory = '8G'

    def end(self):
        super(RmatsStatAgent, self).end()

class RmatsStatTool(Tool):
    def __init__(self, config):
        super(RmatsStatTool, self).__init__(config)
        self.python = 'program/Python/bin/python'
        self.process_rmats_output_py = os.path.join(self.config.PACKAGE_DIR, 'lnc_rna/process_rmats_stat.py')
        self.rmats_output = os.path.join(self.work_dir, 'rmats_output')
        if os.path.isdir(self.rmats_output):
            shutil.rmtree(self.rmats_output)
        shutil.copytree(self.option('result_dir').path, self.rmats_output)

    def run(self):
        super(RmatsStatTool, self).run()
        self.run_process_rmats_stat()
        self.set_output()
        self.end()

    def run_process_rmats_stat(self):
        cmd = '{} {}'.format(self.python, self.process_rmats_output_py)
        cmd += ' -i {}'.format(self.rmats_output)
        cmd += ' -m {}'.format(self.option('significant_diff'))
        cmd += ' -c {}'.format(self.option('significant_value'))
        cmd += ' -p {}'.format(self.option('delta_psi'))
        cmd_name = 'run_process_rmats_stat'
        self.run_code(cmd_name, cmd)

    def run_code(self, cmd_name, cmd, shell=False):
        if shell:
            cmd = self.config.SOFTWARE_DIR + '/' + cmd
        command = self.add_command(cmd_name, cmd, shell=shell)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info('succeed in running {}'.format(cmd_name))
        elif command.return_code is None:
            self.logger.warn('fail to run {}, try again'.format(cmd_name))
            command.rerun()
            self.wait()
            if command.return_code is 0:
                self.logger.info('succeed in rerunning {}'.format(cmd_name))
            else:
                self.set_error('fail to rerun {}, abord'.format(cmd_name))
        else:
            self.set_error('fail to run {}, abord'.format(cmd_name))

    def set_output(self):
        self.logger.info('start set_output at {}'.format(self.__class__.__name__))
        for basename in os.listdir(self.rmats_output):
            source = os.path.join(self.rmats_output, basename)
            link_name = os.path.join(self.output_dir, basename)
            if os.path.exists(link_name):
                os.remove(link_name)
            os.link(source, link_name)
            self.logger.info('succeed in linking {} to {}'.format(source, link_name))
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
            'id': 'rmats_stat_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'tool',
            'name': 'lnc_rna.structure.rmats_stat',
            'instant': False,
            'options': {
                'result_dir': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/lnc_rna/rmats/output/Vit_vs_Con',
                'delta_psi': 0.1,
                'significant_diff': 'pvalue',
                'significant_value': 0.01,
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    unittest.main()
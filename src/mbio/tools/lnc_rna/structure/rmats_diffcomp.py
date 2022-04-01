# -*- coding: utf-8 -*-
# __author__ = 'gudeqing, qinjincheng'

from biocluster.agent import Agent
from biocluster.tool import Tool
import os
import unittest

class RmatsDiffcompAgent(Agent):
    '''
    last_modify: 2019.03.19
    '''
    def __init__(self, parent):
        super(RmatsDiffcompAgent, self).__init__(parent)
        options = [
            {'name': 'rmats_detail_list', 'type': 'infile', 'format': 'lnc_rna.common'},
            {'name': 'significant_value', 'type': 'float', 'default': None},
            {'name': 'delta_psi', 'type': 'float', 'default': None},
            {'name': 'significant_diff', 'type': 'string', 'default': None},
            {'name': 'type_file', 'type': 'infile', 'format': 'lnc_rna.common'},
            {'name': 'diffcomp_txt', 'type': 'outfile', 'format': 'lnc_rna.common'},
        ]
        self.add_option(options)
        self.step.add_steps('rmats_diffcomp')
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def check_options(self):
        self.logger.info('start check_options at {}'.format(self.__class__.__name__))
        for k, v in self._options.items():
            self.logger.debug('{} = {}'.format(k, v.value))
        self.logger.info('finish check_options at {}'.format(self.__class__.__name__))

    def stepstart(self):
        self.step.rmats_diffcomp.start()
        self.step.update()

    def stepfinish(self):
        self.step.rmats_diffcomp.finish()
        self.step.update()

    def set_resource(self):
        self._cpu = 2
        self._memory = '4G'

    def end(self):
        super(RmatsDiffcompAgent, self).end()

class RmatsDiffcompTool(Tool):
    def __init__(self,config):
        super(RmatsDiffcompTool, self).__init__(config)
        self.python = 'program/Python/bin/python'
        self.process_rmats_detail_py = os.path.join(self.config.PACKAGE_DIR, 'lnc_rna/process_rmats_detail.py')
        self.diffcomp_txt = os.path.join(self.work_dir, 'diffcomp.txt')

    def run(self):
        super(RmatsDiffcompTool, self).run()
        self.run_process_rmats_detail()
        self.set_output()
        self.end()

    def run_process_rmats_detail(self):
        cmd = '{} {}'.format(self.python, self.process_rmats_detail_py)
        cmd += ' -i {}'.format(self.option('rmats_detail_list').path)
        cmd += ' -m {}'.format(self.option('significant_diff'))
        cmd += ' -p {}'.format(self.option('delta_psi'))
        cmd += ' -c {}'.format(self.option('significant_value'))
        cmd += ' -t {}'.format(self.option('type_file').path)
        cmd += ' -o {}'.format(self.diffcomp_txt)
        cmd_name = 'run_process_rmats_detail'
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
        source = self.diffcomp_txt
        link_name = os.path.join(self.output_dir, os.path.basename(source))
        if os.path.exists(link_name):
            os.remove(link_name)
        os.link(source, link_name)
        self.option('diffcomp_txt').set_path(link_name)
        self.logger.info('succeed in linking {} to {}'.format(source, link_name))
        self.logger.info('finish set_output at {}'.format(self.__class__.__name__))

class TestFunction(unittest.TestCase):
    '''
    This is test for the tool. Just run this script to do test.
    '''
    def test_fdr(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'rmats_diffcomp_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'tool',
            'name': 'lnc_rna.structure.rmats_diffcomp',
            'instant': False,
            'options': {
                'rmats_detail_list': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/lnc_rna/rmats_diffcomp/rmats_detail.loc2name.list',
                'delta_psi': 0.0,
                'significant_diff': 'fdr',
                'significant_value': 0.05,
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

    def test_pvalue(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'rmats_diffcomp_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'tool',
            'name': 'lnc_rna.structure.rmats_diffcomp',
            'instant': False,
            'options': {
                'rmats_detail_list': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/lnc_rna/rmats_diffcomp/rmats_detail.loc2name.list',
                'significant_value': 0.05,
                'delta_psi': 0.0,
                'significant_diff': 'pvalue',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    unittest.main()